/*
 * Copyright (C) 1993-2017 David Rowe, Brady O'Brien
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package fsk;

import math.Complex;
import math.ComplexMath;
import math.FFT;
import java.util.Arrays;

public final class Demodulator implements Defines {

    private final Complex[] phi_c;
    private final Complex[] samp_old;
    //
    private final double[] fft_est;
    private double norm_rx_timing;
    private double snr;
    //
    private int nin;
    //
    private final FFT fft;
    private final Statistics stats;

    public Demodulator() {
        int i;

        fft = new FFT(FFTSIZE);
        stats = new Statistics();
        samp_old = new Complex[NSTASH];
        nin = FRAME_SAMPLES;   // 32 symbols or 640 samples

        fft_est = new double[FFTSIZE / 2];

        phi_c = new Complex[NUMBER_OF_TONES];

        // Set up rx phases
        for (i = 0; i < NUMBER_OF_TONES; i++) {
            phi_c[i] = ComplexMath.cexp(new Complex());
        }

        for (i = 0; i < NSTASH; i++) {
            samp_old[i] = new Complex();
        }
    }

    public int fsk_get_nin() {
        return nin;
    }

    public double fsk_get_snr() {
        return snr;
    }

    /**
     * Decode the modulation tones over the frame, and provide their FFT bin
     * frequency
     *
     * @param baseband a complex time domain signal of frame samples
     * @param freqs an array of frequencies that have been decoded
     */
    private void frequencyEstimate(Complex[] baseband, double[] freqs) {
        int i, j;

        Complex[] fftdata = new Complex[FFTSIZE];
        double[] magnitude = new double[FFTSIZE / 2];
        double[] result = new double[FFTSIZE / 2];

        Complex dphi = ComplexMath.cexp(new Complex(0.0, TAU / (double) (FFTSIZE - 1)));
        Complex rphi = ComplexMath.timesConjugate(ComplexMath.times(ComplexMath.cexp(new Complex()), 0.5), dphi);

        int f_min = (EST_MIN * FFTSIZE) / SAMPLERATE;
        int f_max = (EST_MAX * FFTSIZE) / SAMPLERATE;
        int f_zero = (EST_SPACE * FFTSIZE) / SAMPLERATE;            // 20

        int fft_samps = FFTSIZE;

        for (j = 0; j < 2; j++) {
            for (i = 0; i < fft_samps; i++) {
                rphi = ComplexMath.times(rphi, dphi);
                // Hann filter
                fftdata[i] = ComplexMath.times(baseband[i + (FFTSIZE * j)], 0.5 - rphi.getReal());
            }

            // on second pass there's only 128 samples
            // so zero out the rest of the FFT input
            if (j == 1) {
                for (i = fft_samps; i < FFTSIZE; i++) {
                    fftdata[i] = new Complex();
                }
            }

            fft.transform(fftdata);

            // Find the magnitude of each freq slot
            for (i = 0; i < (FFTSIZE / 2); i++) {
                magnitude[i] = ComplexMath.absolute(fftdata[i]);
            }

            // Zero out the minimum and maximum ends
            // which will simplify our peak search
            for (i = 0; i < f_min; i++) {
                magnitude[i] = 0.0;
            }

            for (i = (f_max - 1); i < (FFTSIZE / 2); i++) {
                magnitude[i] = 0.0;
            }

            // IIR filter the new frames spectrum with the old
            for (i = 0; i < (FFTSIZE / 2); i++) {
                fft_est[i] = fft_est[i] * (1.0 - TC) + (magnitude[i] * TC);
                result[i] = fft_est[i];
            }

            fft_samps = (nin - FFTSIZE);    // nominally 640 - 512 = 128
        }

        int[] freqi = new int[NUMBER_OF_TONES];
        double max;
        int imax;

        /*
         * At this point we have the Frequency Domain filtered
         * over the 32 symbols (640 PCM samples. So now we can
         * find the peaks and hopefully, if this is a good customer
         * on the other end, we will have the four FSK tones
         */
        for (i = 0; i < NUMBER_OF_TONES; i++) {
            imax = 0;
            max = 0.0;

            for (j = 0; j < (FFTSIZE / 2); j++) {
                if (result[j] >= max) {
                    max = result[j];
                    imax = j;
                }
            }

            // Save the index
            freqi[i] = imax;

            // Blank out Center of tone +/- 20 filters
            f_min = (imax - f_zero);
            f_min = (f_min < 0) ? 0 : f_min;

            f_max = (imax + f_zero);
            f_max = (f_max > (FFTSIZE / 2)) ? (FFTSIZE / 2) : f_max;

            // zero out result for next peak search
            for (j = f_min; j < f_max; j++) {
                result[j] = 0.0;
            }
        }

        // We need the frequencies sorted low to high
        Arrays.sort(freqi);

        // Normalize FFT bins
        for (i = 0; i < NUMBER_OF_TONES; i++) {
            freqs[i] = (double) freqi[i] / (double) FFTSIZE;
        }
    }

    /**
     * Demodulate a complex signal containing 4FSK to a frame of bits
     *
     * @param freqs a utility for the user to show tones detected
     * @param rx_bits an array containing 8 sync + 56 codec bits
     * @param baseband the complex baseband fsk signal
     */
    public void demodulate(double[] freqs, byte[] rx_bits, Complex[] baseband) {
        Complex[][] f_int = new Complex[NUMBER_OF_TONES][(FRAME_SYMBOLS + 1) * OVERSAMPLE];
        Complex[] f_intbuf = new Complex[CYCLES];
        Complex[] phi = new Complex[NUMBER_OF_TONES];
        Complex[] dphi = new Complex[NUMBER_OF_TONES];
        int i, j, m;

        for (i = 0; i < CYCLES; i++) {
            f_intbuf[i] = new Complex();
        }

        for (i = 0; i < NUMBER_OF_TONES; i++) {
            phi[i] = phi_c[i];
        }

        for (i = 0; i < NUMBER_OF_TONES; i++) {
            for (j = 0; j < ((FRAME_SYMBOLS + 1) * OVERSAMPLE); j++) {
                f_int[i][j] = new Complex();
            }
        }

        // Estimate four frequencies from complex frame samples
        frequencyEstimate(baseband, freqs);

        Complex dphi_m;
        int cbuf_i;
        int dc_i;
        int nold = CYCLE_MEMORY - nin;     // 680 - 640 = 40 +/- 10 (nin)
        boolean using_old_samps;

        // Initalize downmixers for each symbol tone
        for (i = 0; i < NUMBER_OF_TONES; i++) {
            // Back the stored phase off to account for re-integraton of old samples
            dphi[i] = ComplexMath.cexp(new Complex(0.0, -TAU * (nold - (CYCLES / OVERSAMPLE)) * freqs[i]));
            phi[i] = ComplexMath.times(phi[i], dphi[i]);
        }

        // Integrate and downsample for symbol tones
        for (m = 0; m < NUMBER_OF_TONES; m++) {
            dphi_m = ComplexMath.cexp(new Complex(0.0, TAU * freqs[m]));
            using_old_samps = true;

            // Pre-fill integration buffer
            for (dc_i = 0; dc_i < (CYCLES - (CYCLES / OVERSAMPLE)); dc_i++) {   // 20 - 2 = 17
                // Switch sample source to new samples when we run out of old ones

                if ((dc_i >= nold) && (using_old_samps == true)) {
                    dc_i = 0;
                    using_old_samps = false;

                    // Recalculate delta-phi after switching to new sample source
                    ComplexMath.normalize(phi[m]);
                }

                if (using_old_samps == true) {
                    f_intbuf[dc_i] = ComplexMath.timesConjugate(samp_old[(NSTASH - nold) + dc_i], phi[m]);
                } else {
                    f_intbuf[dc_i] = ComplexMath.timesConjugate(baseband[dc_i], phi[m]);
                }

                // Spin downconversion phases
                phi[m] = ComplexMath.times(phi[m], dphi_m);
            }

            cbuf_i = dc_i;

            // Integrate over Ts at offsets of Ts/P
            for (i = 0; i < ((FRAME_SYMBOLS + 1) * OVERSAMPLE); i++) {
                // Downconvert and Place Ts/P samples in the integration buffers

                for (j = 0; j < (CYCLES / OVERSAMPLE); j++, dc_i++) {   // 2

                    // Switch sample source to new samples when we run out of old ones
                    if ((dc_i >= nold) && (using_old_samps == true)) {
                        dc_i = 0;
                        using_old_samps = false;

                        // Recalculate delta-phi after switching to new sample source
                        ComplexMath.normalize(phi[m]);
                    }

                    if (using_old_samps == true) {
                        f_intbuf[cbuf_i + j] = ComplexMath.timesConjugate(samp_old[(NSTASH - nold) + dc_i], phi[m]);
                    } else {
                        f_intbuf[cbuf_i + j] = ComplexMath.timesConjugate(baseband[dc_i], phi[m]);
                    }

                    // Spin downconversion phases
                    phi[m] = ComplexMath.times(phi[m], dphi_m);
                }

                cbuf_i += (CYCLES / OVERSAMPLE);    // 2

                if (cbuf_i >= CYCLES) {
                    cbuf_i = 0;
                }

                // Integrate over the integration buffers, save samples
                Complex it = new Complex();

                for (j = 0; j < CYCLES; j++) {
                    it = ComplexMath.add(it, f_intbuf[j]);
                }

                f_int[m][i] = it;
            }
        }

        // Save phases
        for (m = 0; m < NUMBER_OF_TONES; m++) {
            phi_c[m] = phi[m];
        }

        // Stash samples away in the old sample buffer for the next round
        System.arraycopy(baseband, nin - NSTASH, samp_old, 0, NSTASH);

        // Fine Timing Estimation
        // Apply magic nonlinearity to f1_int and f2_int, shift down to 0, extract angle
        // Figure out how much to spin the oscillator to extract magic spectral line
        Complex dphift = ComplexMath.cexp(new Complex(0.0, TAU * ((double) SYMBOLRATE / (double) (OVERSAMPLE * SYMBOLRATE))));
        Complex phi_ft = ComplexMath.cexp(new Complex());
        Complex t_c = new Complex();
        double ft1;

        for (i = 0; i < ((FRAME_SYMBOLS + 1) * OVERSAMPLE); i++) {
            // Get abs^2 of fx_int[i], and add 'em
            ft1 = 0.0;

            for (m = 0; m < NUMBER_OF_TONES; m++) {
                ft1 += ComplexMath.absolute(f_int[m][i]);
            }

            // Down shift and accumulate magic line
            t_c = ComplexMath.add(t_c, ComplexMath.times(phi_ft, ft1));

            // Spin the oscillator for the magic line shift
            phi_ft = ComplexMath.times(phi_ft, dphift);
        }

        // Get the magic angle
        norm_rx_timing = ComplexMath.carg(t_c) / TAU;
        double rx_timing = norm_rx_timing * (double) OVERSAMPLE;

        // Figure out how many samples are needed the next modem cycle
        if (rx_timing > 0.25) {
            nin = FRAME_SAMPLES + (CYCLES / 2);
        } else if (rx_timing < -0.25) {
            nin = FRAME_SAMPLES - (CYCLES / 2);
        } else {
            nin = FRAME_SAMPLES;
        }

        Complex[] t = new Complex[NUMBER_OF_TONES];
        double[] tmax = new double[NUMBER_OF_TONES];
        double max;
        int st;
        int sym;

        // Re-sample the integrators with linear interpolation magic
        int low_sample = (int) Math.floor(rx_timing);
        int high_sample = (int) Math.ceil(rx_timing);

        double fract = rx_timing - (double) low_sample;

        for (i = 0; i < NUMBER_OF_TONES; i++) {
            t[i] = new Complex();
        }

        for (i = 0; i < FRAME_SYMBOLS; i++) {     // 32
            st = (i + 1) * OVERSAMPLE;            // 10 - 160

            for (m = 0; m < NUMBER_OF_TONES; m++) {
                t[m] = ComplexMath.times(f_int[m][st + low_sample], 1.0 - fract);
                t[m] = ComplexMath.add(t[m], ComplexMath.times(f_int[m][st + high_sample], fract));

                // Figure mag^2 of each resampled fx_int
                tmax[m] = ComplexMath.absolute(t[m]);
            }

            max = tmax[0];  // Maximum for figuring correct symbol
            sym = 0;        // Index of maximum

            for (m = 0; m < NUMBER_OF_TONES; m++) {
                if (tmax[m] > max) {
                    max = tmax[m];
                    sym = m;
                }
            }

            // Get the actual bits

            rx_bits[(i * 2)]     = (byte)((sym & 0x2) >>> 1);   // LSB
            rx_bits[(i * 2) + 1] = (byte)(sym & 0x1);           // MSB

            stats.Push(Math.sqrt(max));
        }

        snr = 20.0 * Math.log10(stats.Mean() / stats.StandardDeviation());
    }
}
