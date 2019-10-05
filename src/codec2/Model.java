/*
 * Copyright (C) 1993-2017 David Rowe
 *
 * All rights reserved
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package codec2;

import static codec2.Codec2.N_SAMP;
import static codec2.Codec2.TAU;
import java.util.Random;

import math.Complex;
import math.ComplexMath;

public final class Model {

    private final static double BG_THRESH_DB = 40.0;    // only consider low levels signals for bg_est
    private final static double BG_BETA = 0.1;          // averaging filter constant
    private final static double BG_MARGIN_DB = 6.0;
    private final static int MAX_HARMONIC = 80;         // maximum number of harmonics
    //
    private final double[] A;                           // amplitude of each harmonic
    private final double[] phi;                         // phase of each harmonic
    private double Wo;                                  // fundamental frequency estimate in rad/s
    private int L;                                      // number of harmonics
    private boolean voiced;                             // if this frame is voiced
    //
    private final Random rand;

    /*
     * Class to hold model parameters for one 10 ms frame
     */
    public Model() {
        rand = new Random(System.currentTimeMillis());     // seed the noise generator
        phi = new double[MAX_HARMONIC+1];  // 0..80
        A = new double[MAX_HARMONIC+1];
        Wo = 0.0;
        L = 0;
        voiced = false;
    }

    public void reset() {
        for (int i = 0; i <= MAX_HARMONIC; i++) {
            A[i] = 0.0;
        }
    }

    public boolean getVoiced() {
        return voiced;
    }

    public void setVoiced(boolean val) {
        voiced = val;
    }

    public double getA(int index) {
        return A[index];
    }

    public void setA(int index, double val) {
        A[index] = val;
    }

    public double getPhi(int index) {
        return phi[index];
    }

    public void setPhi(int index, double val) {
        phi[index] = val;
    }

    public int getL() {
        return L;
    }

    public void setL(int val) {
        L = val;
    }

    public double getWo() {
        return Wo;
    }

    public void setWo(double val) {
        Wo = val;
    }

    /*
     * Postfilter to improve sound quality for speech with high levels of
     * background noise. Unlike mixed-excitation models requires no bits to be
     * transmitted to handle background noise.
     */
    private void postfilter(double[] bg_est) {
        double e = 1E-12;

        for (int m = 1; m <= L; m++) {
            e += (A[m] * A[m]);
        }

        e = 10.0 * Math.log10(e / L);

        if ((e < BG_THRESH_DB) && !voiced) {
            bg_est[0] *= (1.0 - BG_BETA) + (e * BG_BETA); // IIR
        }

        double thresh = Math.pow(10.0, (bg_est[0] + BG_MARGIN_DB) / 20.0);

        if (voiced == true) {
            for (int m = 1; m <= L; m++) {
                if (A[m] < thresh) {
                    phi[m] = TAU * rand.nextDouble(); // random value 0.0..TAU
                }
            }
        }
    }

    /*
     * Synthesizes phases based on SNR and a rule based approach. No phase
     * parameters are required apart from the SNR (which can be reduced to a 1
     * bit V/UV decision per frame).
     */
    public void phase_synth_zero_order(double[] ex_phase, double[] bg_est, Complex[] Aw) {
        Complex Ex;

        /*
         * Update excitation fundamental phase track, this sets the position
         * of each pitch pulse during voiced speech.
         */
        ex_phase[0] += (Wo * N_SAMP);
        ex_phase[0] -= TAU * Math.floor(ex_phase[0] / TAU + 0.5);

        for (int m = 1; m <= L; m++) {

            /*
             * generate excitation
             */
            if (voiced == true) {
                Ex = ComplexMath.cexp(new Complex(0.0, (ex_phase[0] * (double) m)));
            } else {
                /*
                 * When a few samples were tested I found that LPC filter
                 * phase is not needed in the un-voiced case, but no harm in
                 * keeping it.
                 */

                Ex = ComplexMath.cexp(new Complex(0.0, TAU * rand.nextDouble()));
            }

            /*
             * filter using LPC filter
             * modify sinusoidal phase
             */
            phi[m] = ComplexMath.carg(ComplexMath.times(Aw[m], Ex));
        }
        
        postfilter(bg_est);
    }
}
