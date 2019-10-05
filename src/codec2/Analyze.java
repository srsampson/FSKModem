/*
 * Copyright (C) 1993-2017 David Rowe
 *
 * All rights reserved
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package codec2;

import static codec2.Codec2.FFT_SIZE;
import static codec2.Codec2.FS;
import static codec2.Codec2.M_PITCH;
import static codec2.Codec2.NW;
import static codec2.Codec2.N_SAMP;
import static codec2.Codec2.P_MAX;
import static codec2.Codec2.P_MIN;
import static codec2.Codec2.TAU;

import math.Complex;
import math.ComplexMath;
import math.FFT;

public final class Analyze {

    private static final double VOICING_THRESHOLD_DB = 6.0;
    private static final double R = TAU / FFT_SIZE;
    private static final double WO_MAX = TAU / P_MIN;
    private static final double WO_MIN = TAU / P_MAX;
    private static final double CNLP = 0.3;         // post processor constant
    private static final double COEFF = 0.95;       // notch filter parameter
    //
    private static final int NLP_NTAP = 48;         // Decimation LPF order
    private static final int DEC = 5;               // decimation factor
    //
    private static final int MIN_BIN = FFT_SIZE * DEC / P_MAX;   // 16
    private static final double PREV_BIN = (4000.0 / Math.PI) * (double)(FFT_SIZE * DEC) / FS;
    //
    private final FFT fft;
    //
    private final Complex[] fw;             // DFT of squared signal (input)
    //
    private final double[] W;
    private final double[] mem_fir;         // decimation FIR filter memory
    private final double[] sq;              // squared speech samples
    private final double[] coswintab;       // optimization
    private final double[] Sn;
    //
    private final double[] hamming;
    //
    private double mem_x;
    private double mem_y;
    //
    private double previous_Wo;

    /*
     * 48 tap 600Hz low pass FIR filter coefficients
     */
    private static final double[] NLPCOEF = {
        -1.0818124e-03,
        -1.1008344e-03,
        -9.2768838e-04,
        -4.2289438e-04,
        5.5034190e-04,
        2.0029849e-03,
        3.7058509e-03,
        5.1449415e-03,
        5.5924666e-03,
        4.3036754e-03,
        8.0284511e-04,
        -4.8204610e-03,
        -1.1705810e-02,
        -1.8199275e-02,
        -2.2065282e-02,
        -2.0920610e-02,
        -1.2808831e-02,
        3.2204775e-03,
        2.6683811e-02,
        5.5520624e-02,
        8.6305944e-02,
        1.1480192e-01,
        1.3674206e-01,
        1.4867556e-01,
        //
        1.4867556e-01,
        1.3674206e-01,
        1.1480192e-01,
        8.6305944e-02,
        5.5520624e-02,
        2.6683811e-02,
        3.2204775e-03,
        -1.2808831e-02,
        -2.0920610e-02,
        -2.2065282e-02,
        -1.8199275e-02,
        -1.1705810e-02,
        -4.8204610e-03,
        8.0284511e-04,
        4.3036754e-03,
        5.5924666e-03,
        5.1449415e-03,
        3.7058509e-03,
        2.0029849e-03,
        5.5034190e-04f,
        -4.2289438e-04,
        -9.2768838e-04,
        -1.1008344e-03,
        -1.0818124e-03
    };

    public Analyze(FFT fftencode) {
        int i;

        fft = fftencode;
        hamming = new double[M_PITCH];
        W = new double[FFT_SIZE];      // changed to double, as imag is always 0.0
        new Hamming(fftencode, hamming, W);

        mem_fir = new double[NLP_NTAP];
        coswintab = new double[M_PITCH / DEC];
        sq = new double[M_PITCH];
        Sn = new double[M_PITCH];
        previous_Wo = 0.0;

        for (i = 0; i < M_PITCH; i++) {
            Sn[i] = 1.0;
        }

        mem_x = 0.0;
        mem_y = 0.0;

        fw = new Complex[FFT_SIZE];

        for (i = 0; i < FFT_SIZE; i++) {
            fw[i] = new Complex();
        }

        /*
         * Pre-compute decimation window
         */
        for (i = 0; i < (M_PITCH / DEC); i++) {
            coswintab[i] = 0.5 - 0.5 * Math.cos(TAU * i / (double)(M_PITCH / DEC - 1));
        }
    }

    public void reset() {
        for (int i = 0; i < M_PITCH; i++) {
            Sn[i] = 1.0;
        }
    }

    /*
     * Determines the pitch in samples using a Non Linear Pitch (NLP) algorithm.
     *
     * Returns the pitch in Hz.
     */
    private double nlp() {
        double notch;
        double thresh, lmax;
        int b, bmin, bmax, lmax_bin;
        int i, j;

        /*
         * Square, notch filter at DC, and LP filter vector
         */
        for (i = (M_PITCH - N_SAMP); i < M_PITCH; i++) {
            sq[i] = (Sn[i] * Sn[i]);

            notch = sq[i] - mem_x;
            notch += COEFF * mem_y;

            mem_x = sq[i];
            mem_y = notch;
            sq[i] = notch + 1.0;

            System.arraycopy(mem_fir, 1, mem_fir, 0, (NLP_NTAP - 1));   // shift memory left
            mem_fir[NLP_NTAP - 1] = sq[i];

            sq[i] = 0.0;
            for (j = 0; j < NLP_NTAP; j++) {
                sq[i] += mem_fir[j] * NLPCOEF[j];
            }
        }

        // decimate and DFT
        for (i = 0; i < (M_PITCH / DEC); i++) {
            fw[i] = new Complex(sq[i * DEC] * coswintab[i], 0.0);
        }

        // fill the rest with 0.0
        for (i = (M_PITCH / DEC); i < FFT_SIZE; i++) {
            fw[i] = new Complex();
        }

        fft.transform(fw);

        for (i = 0; i < FFT_SIZE; i++) {
            fw[i] = new Complex(ComplexMath.square(fw[i]), 0.0);
        }

        /* find global peak */
        double gmax = 0.0;
        int gmax_bin = FFT_SIZE * DEC / P_MAX;      // 16

        /*
         * 512 * 5 / 160(P_MAX) = 16  and 512 * 5 / 20(P_MIN) = 128
         */
        for (i = FFT_SIZE * DEC / P_MAX; i <= FFT_SIZE * DEC / P_MIN; i++) {
            if (fw[i].getReal() > gmax) {
                gmax = fw[i].getReal();
                gmax_bin = i;
            }
        }

        /* Shift samples in buffer to make room for new samples */
        System.arraycopy(sq, N_SAMP, sq, 0, M_PITCH - N_SAMP);

        int prev_fo_bin = (int) (previous_Wo * PREV_BIN);

        /*
         * post process estimate by searching submultiples
         */
        int mult = 2;
        int cmax_bin = gmax_bin;

        while ((gmax_bin / mult) >= MIN_BIN) {
            b = (gmax_bin / mult);			// determine search interval
            bmin = (int) (0.8 * b);
            bmax = (int) (1.2 * b);

            if (bmin < MIN_BIN) {
                bmin = MIN_BIN;
            }

            /*
             * lower threshold to favour previous frames pitch estimate,
             * this is a form of pitch tracking
             */
            if ((prev_fo_bin > bmin) && (prev_fo_bin < bmax)) {
                thresh = CNLP * 0.5 * gmax;
            } else {
                thresh = CNLP * gmax;
            }

            lmax = 0.0;
            lmax_bin = bmin;

            for (b = bmin; b <= bmax; b++) {        // look for maximum in interval
                if (fw[b].getReal() > lmax) {
                    lmax = fw[b].getReal();
                    lmax_bin = b;
                }
            }

            if (lmax > thresh) {
                if ((lmax > fw[lmax_bin - 1].getReal()) && (lmax > fw[lmax_bin + 1].getReal())) {
                    cmax_bin = lmax_bin;
                }
            }

            mult++;
        }

        // pitch = sample rate / best Fo
        return FS / ((double) cmax_bin * FS / (FFT_SIZE * DEC));
    }

    /*
     * Extract sinusoidal model parameters from input speech samples. 10
     * milliseconds of new speech samples are added during each call.
     */
    public void analyze_one_frame(Model model, short[] speech, int index) {
        Complex[] Swe = new Complex[FFT_SIZE];   // complex representation
        int i;

        /*
         * Initialize the working buffer with complex zero's
         */
        for (i = 0; i < FFT_SIZE; i++) {
            Swe[i] = new Complex();
        }

        /*
         * Prepare for new 80 samples.
         *
         * The Sn array is initialized to all 1.0 values in reset()
         */
        System.arraycopy(Sn, N_SAMP, Sn, 0, M_PITCH - N_SAMP);    // M = 320, N = 80, M-N = 240

        /*
         * Now add the new samples to the end
         */
        for (i = 0; i < N_SAMP; i++) {                       // N = 80
            Sn[i + (M_PITCH - N_SAMP)] = (double) speech[i + index];
        }

        /*
         * Center analysis window on time axis. We need to arrange input
         * time domain to make FFT phases correct
         */
        // move 2nd half to start of FFT input vector (160,299 --> 0,139)
        for (i = 0; i < NW / 2; i++) {  // NW/2 = 140
            Swe[i] = new Complex(Sn[i + M_PITCH / 2] * hamming[i + M_PITCH / 2], 0.0);
        }

        // move 1st half to end of FFT input vector (20,159 --> 372,511) 
        for (i = 0; i < NW / 2; i++) {
            Swe[(FFT_SIZE - NW / 2 + i)] = new Complex(Sn[i + M_PITCH / 2 - NW / 2]
                    * hamming[i + M_PITCH / 2 - NW / 2], 0.0);
        }

        fft.transform(Swe);       // now in frequency domain

        double Wo = TAU / nlp();

        model.setWo(Wo);
        model.setL((int) (Math.PI / Wo));


        // calculate model parameters
        two_stage_pitch_refinement(model, Swe);   // refine the pitch
        estimate_amplitudes(model, Swe);
        est_voicing_mbe(model, Swe);

        previous_Wo = model.getWo();         // used in pitch on next pass
    }

    /*
     * Refines the current pitch estimate using the harmonic sum pitch
     * estimation technique.
     */
    private void two_stage_pitch_refinement(Model model, Complex[] Swe) {
        double pmin, pmax;	// pitch refinement minimum, maximum

        // Coarse refinement
        pmax = TAU / model.getWo() + 5.0;
        pmin = TAU / model.getWo() - 5.0;
        hs_pitch_refinement(model, Swe, pmin, pmax, 1.0);    // step 1.0 course

        // Fine refinement
        pmax = TAU / model.getWo() + 1.0;
        pmin = TAU / model.getWo() - 1.0;
        hs_pitch_refinement(model, Swe, pmin, pmax, 0.25);   // step 0.25 fine

        double Wo = model.getWo();

        // Limit range
        if (Wo < WO_MIN) {
            Wo = WO_MIN;
            model.setWo(Wo);
        } else if (Wo > WO_MAX) {    // changed to if-else
            Wo = WO_MAX;
            model.setWo(Wo);
        }

        model.setL((int) Math.floor(Math.PI / Wo));   // was floor()
    }

    /*
     * Harmonic sum pitch refinement function.
     *
     * pmin pitch search range minimum pmax pitch search range maximum step
     * pitch search step size model current pitch estimate in model.Wo
     *
     * model refined pitch estimate in model.Wo
     */
    private void hs_pitch_refinement(Model model, Complex[] Swe, double pmin, double pmax, double pstep) {
        double E;	// energy for current pitch
        double Wo;	// current "test" fundamental freq.
        double Wom;	// Wo that maximizes E
        double Em;	// maximum energy
        double p;	// current pitch
        int m;		// loop variable
        int b;		// bin for current harmonic center

        Wom = model.getWo();
        model.setL((int) (Math.PI / Wom));	// use initial pitch est. for L

        Em = 0.0;

        // Determine harmonic sum for a range of Wo values
        for (p = pmin; p <= pmax; p += pstep) {
            E = 0.0;
            Wo = TAU / p;

            // Sum harmonic magnitudes
            for (m = 1; m <= model.getL(); m++) {
                b = (int) (m * Wo / R + 0.5);
                E += ComplexMath.square(Swe[b]);
            }

            // Compare to see if this is a maximum
            if (E > Em) {
                Em = E;
                Wom = Wo;
            }
        }

        model.setWo(Wom);
    }

    /*
     * Estimates the complex amplitudes of the harmonics.
     */
    private void estimate_amplitudes(Model model, Complex[] Swe) {
        double den;
        int am, bm, b, i, m;

        double tmp = model.getWo() / R;
        int l = model.getL();

        for (m = 1; m <= l; m++) {
            am = (int) ((m - 0.5) * tmp + 0.5);     // lower
            bm = (int) ((m + 0.5) * tmp + 0.5);     // upper
            b = (int) (m * tmp + 0.5);              // center

            // Estimate ampltude of harmonic
            den = 0.0;

            for (i = am; i < bm; i++) {
                den += ComplexMath.square(Swe[i]);
            }

            model.setA(m, Math.sqrt(den));
        }
    }

    /*
     * Returns the error of the MBE cost function for a given F0.
     *
     * Note: I think a lot of the operations below can be simplified as
     * W[].getImaginary() = 0 and has been normalized such that den always
     * equals 1.
     */
    private void est_voicing_mbe(Model model, Complex[] Swe) {
        Complex Sw;
        Complex Am;
        Complex Ew;
        double den;
        int i, m, al, bl, offset;

        double signal = 1E-4;
        int l = model.getL();

        for (i = 1; i <= (l / 4); i++) {
            signal += (model.getA(i) * model.getA(i));
        }

        double Wo = model.getWo();
        double error = 1E-4;

        /* Just test across the harmonics in the first 1000 Hz (L/4) */
        for (i = 1; i <= (l / 4); i++) {
            Am = new Complex();
            den = 0.0;

            al = (int) Math.ceil((i - 0.5) * Wo * FFT_SIZE / TAU);
            bl = (int) Math.ceil((i + 0.5) * Wo * FFT_SIZE / TAU);

            offset = (int) (FFT_SIZE / 2.0 - i * Wo * FFT_SIZE / TAU + 0.5);

            for (m = al; m < bl; m++) {
                Am = ComplexMath.add(Am, ComplexMath.times(Swe[m], W[offset + m]));
                den += W[offset + m] * W[offset + m];
            }

            Am = ComplexMath.divide(Am, den);

            /*
             * Determine error between estimated harmonic and original
             */
            for (m = al; m < bl; m++) {
                Sw = ComplexMath.times(Am, W[offset + m]);
                Ew = new Complex(Swe[m].getReal() - Sw.getReal(),
                        Swe[m].getImaginary() - Sw.getImaginary());
                error += ComplexMath.square(Ew);
            }
        }

        if (10.0 * Math.log10(signal / error) > VOICING_THRESHOLD_DB) {
            model.setVoiced(true);
        } else {
            model.setVoiced(false);
        }

        /*
         * post processing, helps clean up some voicing errors
         *
         * Determine the ratio of low freqency to high frequency energy,
         * voiced speech tends to be dominated by low frequency energy,
         * unvoiced by high frequency. This measure can be used to
         * determine if we have made any gross errors.
         */
        double elow = 1E-4;
        double ehigh = 1E-4;

        for (i = 1; i <= (l / 2); i++) {
            elow += (model.getA(i) * model.getA(i));
        }

        for (i = (l / 2); i <= l; i++) {
            ehigh += (model.getA(i) * model.getA(i));
        }

        double eratio_dB = 10.0 * Math.log10(elow / ehigh);

        if (model.getVoiced() == false) {
            /*
             * Look for Type 1 errors, strongly Voiced speech that has been
             * accidentally declared Un-Voiced
             */

            if (eratio_dB > 10.0) {
                model.setVoiced(true);
            }
        } else if (model.getVoiced() == true) {
            /*
             * Look for Type 2 errors, strongly Un-Voiced speech that has been
             * accidentally declared Voiced
             */

            if (eratio_dB < -10.0) {
                model.setVoiced(false);
            }

            /*
             * A common source of Type 2 errors is the pitch estimator
             * gives a low (50Hz) estimate for Un-Voiced speech, which gives a
             * good match with noise due to the close harmoonic spacing.
             * These errors are much more common than people with 50Hz
             * pitch, so we have just a small eratio_dB threshold.
             */
            if ((eratio_dB < -4.0) && (model.getWo() <= ((60.0 * TAU) / FS))) {
                model.setVoiced(false);
            }
        }
    }
}
