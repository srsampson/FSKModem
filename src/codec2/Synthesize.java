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
import static codec2.Codec2.N_SAMP;
import static codec2.Codec2.TAU;

public final class Synthesize {

    private final FFT fft;
    //
    private final Complex[] Swd;        // the decode speech spectrum
    //
    private final double[] Sn_;
    private final double[] bg_est;      // background noise estimate for post filter
    private final double[] ex_phase;    // phase tracking
    private final double[] Pn;

    public Synthesize(FFT fftdecode) {
        fft = fftdecode;
        Swd = new Complex[FFT_SIZE];
        Sn_ = new double[N_SAMP * 2];
        Pn = new double[N_SAMP * 2];
        bg_est = new double[1];
        ex_phase = new double[1];

        new Parzen(Pn);
    }

    public void reset() {
        bg_est[0] = 0.0;
        ex_phase[0] = 0.0;
    }

    /*
     * Synthesize 80 speech samples (10ms) from model parameters. Limits output
     * level to protect ears when there are bit errors or the input is over
     * driven.
     *
     * This doesn't correct or mask bit errors, just reduces the effects.
     */
    public void synthesize_one_frame(Model model, short[] speech, int index, Complex[] Aw) {
        double gain;
        int i;

        model.phase_synth_zero_order(ex_phase, bg_est, Aw);
        synthesize(model);                  // get updated Sn_

        /*
         * find maximum sample in frame for ear protection
         */
        double max_sample = 0.0;

        for (i = 0; i < N_SAMP; i++) {
            if (Sn_[i] > max_sample) {
                max_sample = Sn_[i];
            }
        }

        /*
         * determine how far above set point
         */
        double over = max_sample / 30000.0;

        /*
         * If we are x dB over set point we reduce level by 2x dB, this
         * attenuates major excursions in amplitude (likely to be caused
         * by bit errors) more than smaller ones
         */
        if (over > 1.0) {
            gain = (over * over);

            for (i = 0; i < N_SAMP; i++) {
                Sn_[i] /= gain;
            }
        }

        for (i = 0; i < N_SAMP; i++) {           // 80
            if (Sn_[i] > 32767.0) {
                speech[i + index] = 32760;
            } else if (Sn_[i] < -32767.0) {
                speech[i + index] = -32760;
            } else {
                speech[i + index] = (short) Sn_[i];
            }
        }
    }

    /*
     * Synthesize a speech signal in the frequency domain from the sinusoidal
     * model parameters. Uses overlap-add with a trapezoidal window to smoothly
     * interpolate between frames.
     */
    private void synthesize(Model model) {
        int i, j;

        /*
         * Initialize the working buffer to complex zero
         */
        for (i = 0; i < FFT_SIZE; i++) {
            Swd[i] = new Complex();
        }

        // Update window by shifting 80 samples left (10 ms)
        System.arraycopy(Sn_, N_SAMP, Sn_, 0, N_SAMP - 1);  // 0 - 78
        Sn_[N_SAMP - 1] = 0.0;                              // 79

        /*
         * Now set up frequency domain synthesized speech
         */
        double tmp = model.getWo() * FFT_SIZE / TAU;

        for (i = 1; i <= model.getL(); i++) {
            int b = (int) (i * tmp + 0.5);

            if (b > ((FFT_SIZE / 2) - 1)) {
                b = (FFT_SIZE / 2) - 1;
            }

            Swd[b] = ComplexMath.times(ComplexMath.cexp(new Complex(0.0, model.getPhi(i))), model.getA(i));
            Swd[FFT_SIZE - b] = ComplexMath.conjugate(Swd[b]);
        }

        // Perform inverse DFT
        fft.itransform(Swd);

        // Overlap add to previous samples
        for (i = 0; i < (N_SAMP - 1); i++) {
            Sn_[i] += (Swd[FFT_SIZE - N_SAMP + 1 + i].getReal() * Pn[i]);
        }

        // put the new data on the end of the window
        for (i = N_SAMP - 1, j = 0; i < (N_SAMP * 2); i++, j++) {
            Sn_[i] = (Swd[j].getReal() * Pn[i]);
        }
    }
}
