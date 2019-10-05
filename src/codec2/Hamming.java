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
import static codec2.Codec2.M_PITCH;
import static codec2.Codec2.NW;
import static codec2.Codec2.TAU;

import math.Complex;
import math.FFT;

/*
 * This whole class should be garbage collected after run,
 * since it will never be used again
 */
public final class Hamming {

    private final FFT fft;
    //
    private Complex[] W;
    
    public Hamming(FFT fftencode, double[] hamming, double[] w) {
        int i, j;

        fft = fftencode;

        /*
         * Generate Hamming window centered on M-sample pitch analysis window
         *
         *         0            M/2           M-1
         *         |-------------|-------------|
         *               |-------|-------|
         *                  NW samples
         *
         * All our analysis/synthsis is centred on the M/2 sample.
         */
        for (i = 0; i < (M_PITCH / 2) - (NW / 2); i++) {
            hamming[i] = 0.0;
        }

        double m = 0.0;

        for (i = (M_PITCH / 2) - (NW / 2), j = 0; i < (M_PITCH / 2) + (NW / 2); i++, j++) {
            hamming[i] = 0.5 - 0.5 * Math.cos(TAU * j / (NW - 1));     // actually a hanning
            m += (hamming[i] * hamming[i]);
        }

        for (i = (M_PITCH / 2) + (NW / 2); i < M_PITCH; i++) {
            hamming[i] = 0.0;
        }

        /*
         * Normalize - makes freq domain amplitude estimation straight forward
         */
        m = Math.sqrt(m * FFT_SIZE);

        for (i = 0; i < M_PITCH; i++) {
            hamming[i] /= m;
        }

        /*
         * Generate DFT of analysis window, used for later processing.  Note
         * we modulo FFT_SIZE shift the time domain window w[], this makes the
         * imaginary part of the DFT W[] equal to zero as the shifted w[] is
         * even about the n=0 time axis if NW is odd.  Having the imag part
         * of the DFT W[] makes computation easier.
         *
         *          0                      FFT_SIZE-1
         *          |-------------------------|
         *
         *           ----\               /----
         *                \             /
         *                 \           /      <- shifted version of window w[n]
         *                  \         /
         *                   \       /
         *                    -------
         *
         *          |---------|     |---------|
         *            NW/2              NW/2
         */
        
        W = new Complex[FFT_SIZE];

        for (i = 0; i < FFT_SIZE; i++) {
            W[i] = new Complex();
        }

        for (i = 0; i < NW / 2; i++) {
            W[i] = new Complex(hamming[i + M_PITCH / 2], 0.0);
        }

        for (i = FFT_SIZE - NW / 2, j = M_PITCH / 2 - NW / 2; i < FFT_SIZE; i++, j++) {
            W[i] = new Complex(hamming[j], 0.0);
        }

        fft.transform(W);

        /*
         * Re-arrange W[] to be symmetrical about FFT_SIZE/2.  Makes later
         * analysis convenient.
         *
         * Before:
         *
         *
         *          0                 FFT_SIZE-1
         *          |----------|---------|
         *          __                   _
         *            \                 /
         *             \_______________/
         *
         * After:
         *
         *          0                 FFT_SIZE-1
         *          |----------|---------|
         *                    ___
         *                   /   \
         *          ________/     \_______
         *
         */

        for (i = 0; i < FFT_SIZE / 2; i++) {
            w[i] = W[i + FFT_SIZE / 2].getReal();
            w[i + FFT_SIZE / 2] = W[i].getReal();
        }
    }
}
