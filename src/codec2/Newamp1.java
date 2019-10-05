/*
 * Copyright (C) 1993-2017 David Rowe
 *
 * All rights reserved
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package codec2;

import static codec2.Codec2.FS;
import static codec2.Codec2.MAX_AMP;
import static codec2.Codec2.NEWAMP1_K;
import static codec2.Codec2.NEWAMP1_PHASE_NFFT;
import static codec2.Codec2.NS;
import static codec2.Codec2.P_MAX;
import static codec2.Codec2.P_MIN;
import static codec2.Codec2.TAU;

import math.Complex;
import math.ComplexMath;
import math.FFT;

public final class Newamp1 {

    private static final int MBEST_ENTRIES = 5;
    private static final int MBEST_STAGES = 4;
    //
    private static final double WO_MIN = Math.log10((TAU / P_MAX));
    private static final double WO_MAX = Math.log10((TAU / P_MIN));
    //
    private static final double MEL200 = Math.floor(2595.0 * Math.log10(1.0 + 200 / 700.0) + 0.5);
    private static final double MEL3700 = Math.floor(2595.0 * Math.log10(1.0 + 3700 / 700.0) + 0.5);
    private static final double MELSTEP = (MEL3700 - MEL200);
    private static final double PHASE_SCALE = (20.0 / Math.log(10.0));
    //
    private final Codebooknewamp1 codebook_newamp;
    private final Codebooknewamp1_energy codebook_energy;
    private final FFT fft;
    //
    private final double[] rate_K_sample_freqs_kHz;
    private final double[] rate_K_vec;

    private class MBest {

        private final int[] index;
        private double error;

        public MBest() {
            index = new int[MBEST_STAGES];
            error = 1E32;
        }

        public void reset() {
            for (int i = 0; i < MBEST_STAGES; i++) {
                index[i] = 0;
            }

            error = 1E32;
        }

        public void setIndex(int val, int ind) {
            index[ind] = val;
        }

        public int getIndex(int ind) {
            return index[ind];
        }

        public void setError(double val) {
            error = val;
        }

        public double getError() {
            return error;
        }
    }

    public Newamp1(FFT fftphase) {
        fft = fftphase;
        rate_K_sample_freqs_kHz = new double[NEWAMP1_K];
        rate_K_vec = new double[NEWAMP1_K];

        double step = MELSTEP / (double) (NEWAMP1_K - 1);
        double mel = MEL200;

        codebook_newamp = new Codebooknewamp1();
        codebook_energy = new Codebooknewamp1_energy();

        for (int k = 0; k < NEWAMP1_K; k++) {
            rate_K_sample_freqs_kHz[k] = 0.7 * (Math.pow(10.0, (mel / 2595.0)) - 1.0);
            mel += step;
        }
    }

    public void newamp1_indexes_to_model(Model[] model_, double[][] interpolated_surface_,
            double[] Wo_left, boolean[] voicing_left, double[] prev_rate_K_vec_,
            Complex[][] H, int[] indexes) {
        double[] rate_K_vec_ = new double[NEWAMP1_K];
        double[] left_vec;
        double[] right_vec;
        int i, k;
        double[] aWo_ = new double[4];
        boolean[] avoicing_ = new boolean[4];
        int[] aL_ = new int[4];
        double Wo_right;
        boolean voicing_right;

        // extract latest rate K vector
        newamp1_indexes_to_rate_K_vec(rate_K_vec_, indexes);

        // decode latest Wo and voicing
        if (indexes[3] != 0) {
            Wo_right = decode_log_Wo(indexes[3], 6);
            voicing_right = true;
        } else {
            Wo_right = TAU / 100.0;
            voicing_right = false;
        }

        left_vec = prev_rate_K_vec_;
        right_vec = rate_K_vec_;

        // interpolate 25Hz rate K vec back to 100Hz
        newamp1_interpolate(interpolated_surface_, left_vec, right_vec);

        /* interpolate 25Hz v and Wo back to 100Hz */
        interp_Wo_v(aWo_, aL_, avoicing_, Wo_left[0], Wo_right, voicing_left[0], voicing_right);

        // back to rate L amplitudes, synthesis phase for each frame
        for (i = 0; i < 4; i++) {
            model_[i].setWo(aWo_[i]);
            model_[i].setL(aL_[i]);
            model_[i].setVoiced(avoicing_[i]);

            resample_rate_L(model_[i], interpolated_surface_, i);
            determine_phase(H, i, model_[i]);
        }

        // update memories for next time
        for (k = 0; k < NEWAMP1_K; k++) {
            prev_rate_K_vec_[k] = rate_K_vec_[k];
        }

        Wo_left[0] = Wo_right;
        voicing_left[0] = voicing_right;
    }

    public void newamp1_model_to_indexes(Model model, int[] indexes) {
        double[] rate_K_vec_no_mean = new double[NEWAMP1_K];
        double[] rate_K_vec_no_mean_ = new double[NEWAMP1_K];
        double[] mean = new double[1];
        int k;

        // convert variable rate L to fixed rate K
        resample_const_rate_f(model, rate_K_vec);

        // remove mean and two stage VQ
        double sum = 0.0;

        for (k = 0; k < NEWAMP1_K; k++) {
            sum += rate_K_vec[k];
        }

        mean[0] = sum / NEWAMP1_K;

        for (k = 0; k < NEWAMP1_K; k++) {
            rate_K_vec_no_mean[k] = rate_K_vec[k] - mean[0];
        }

        rate_K_mbest_encode(indexes, rate_K_vec_no_mean, rate_K_vec_no_mean_);

        // scalar quantise mean (effectively the frame energy)
        indexes[2] = quantize(codebook_energy.getCodebook(0).getCodeBookArray(), mean, codebook_energy.getCodebook(0).getNumberOfElements());

        // scalar quantise Wo.  We steal the smallest Wo index
        // to signal an unvoiced frame
        if (model.getVoiced()) {
            int index = encode_log_Wo(model.getWo(), 6);

            if (index == 0) {
                index = 1;
            }

            indexes[3] = index;
        } else {
            indexes[3] = 0;
        }
    }

    /*
     * Insert the results in a vector for codebook entry comparison. The
     * list is ordered by error, so those entries with the smallest error
     * will be first on the list.
     */
    private void mbest_insert(MBest[] mbest, int[] index, double error) {
        int i, j;

        for (i = 0; i < MBEST_ENTRIES; i++) {
            if (error < mbest[i].getError()) {
                for (j = MBEST_ENTRIES - 1; j > i; j--) {
                    mbest[j] = mbest[j - 1];
                }

                for (j = 0; j < MBEST_STAGES; j++) {
                    mbest[i].setIndex(index[j], j);
                }

                mbest[i].setError(error);
                break;
            }
        }
    }

    /*
     * Searches vec[] to a codebook of vectors, and maintains a list of the mbest
     * closest matches.
     */
    private void mbest_search(float[] cb, double[] vec, int m, MBest[] mbest, int[] index) {
        double e, diff;
        int i, j;

        for (j = 0; j < m; j++) {
            e = 0.0;

            for (i = 0; i < NEWAMP1_K; i++) {
                diff = (double) cb[j * NEWAMP1_K + i] - vec[i];
                e += (diff * diff);
            }

            index[0] = j;
            mbest_insert(mbest, index, e);
        }
    }

    private void mag_to_phase(double[] phase, double[] Gdbfk) {
        Complex[] Sdb = new Complex[NEWAMP1_PHASE_NFFT];
        Complex[] cf = new Complex[NEWAMP1_PHASE_NFFT];
        int i;

        // stupid java tricks
        for (i = 0; i < NEWAMP1_PHASE_NFFT; i++) {
            Sdb[i] = new Complex();
            cf[i] = new Complex();
        }

        /* install negative frequency components, 1/Nfft takes into
            account kiss fft lack of scaling on ifft */
        Sdb[0] = new Complex(Gdbfk[0], 0.0);

        for (i = 1; i < NS; i++) {
            Sdb[i] = new Complex(Gdbfk[i], 0.0);
            Sdb[NEWAMP1_PHASE_NFFT - i] = new Complex(Gdbfk[i], 0.0);
        }

        // compute real cepstrum from log magnitude spectrum
        fft.itransform(Sdb);

        for (i = 0; i < NEWAMP1_PHASE_NFFT; i++) {
            Sdb[i] = ComplexMath.divide(Sdb[i], (double) NEWAMP1_PHASE_NFFT);
        }

        // Fold cepstrum to reflect non-min-phase zeros inside unit circle
        cf[0] = Sdb[0];

        for (i = 1; i < NS - 1; i++) {
            cf[i] = ComplexMath.add(Sdb[i], Sdb[NEWAMP1_PHASE_NFFT - i]);
        }

        cf[NS - 1] = Sdb[NS - 1];

        // cf = dB_magnitude + j * minimum_phase
        fft.transform(cf);

        /*
         * The maths says we are meant to be using log(x), not 20*log10(x), so
         * we need to scale the phase to account for this: log(x) =
         * 20*log10(x)/scale
         */
        for (i = 0; i < NS; i++) {
            phase[i] = cf[i].getImaginary() / PHASE_SCALE;
        }
    }

    /*
     * We add some variables here to translate the C pointers into indexes.
     */
    private void interp_para(double[] y, int yindex, double[] xp, int xpindex,
            double[] yp, int ypindex, int np, double[] x, int xindex, int n) {
        double xi, x1, y1, x2, y2, x3, y3, a, b;
        int k = 0;

        for (int i = 0; i < n; i++) {
            xi = x[i + xindex];

            /* k is index into xp of where we start 3 points used to form parabola */
            while ((xp[k + xpindex + 1] < xi) && (k < (np - 3))) {
                k++;
            }

            x1 = xp[k + xpindex];
            y1 = yp[k + ypindex];

            x2 = xp[k + xpindex + 1];
            y2 = yp[k + ypindex + 1];

            x3 = xp[k + xpindex + 2];
            y3 = yp[k + ypindex + 2];

            a = ((y3 - y2) / (x3 - x2) - (y2 - y1) / (x2 - x1)) / (x3 - x1);
            b = ((y3 - y2) / (x3 - x2) * (x2 - x1) + (y2 - y1) / (x2 - x1) * (x3 - x2)) / (x3 - x1);

            y[i + yindex] = a * (xi - x2) * (xi - x2) + b * (xi - x2) + y2;
        }
    }

    private void resample_const_rate_f(Model model, double[] rate_K_vec) {
        double[] AmdB = new double[MAX_AMP + 1];
        double[] rate_L_sample_freqs_kHz = new double[MAX_AMP + 1];
        double AmdB_peak;
        int m;

        /* convert rate L=pi/Wo amplitude samples to fixed rate K */
        AmdB_peak = -100.0;

        for (m = 1; m <= model.getL(); m++) {
            AmdB[m] = 20.0 * Math.log10(model.getA(m) + 1E-16);

            if (AmdB[m] > AmdB_peak) {
                AmdB_peak = AmdB[m];
            }

            rate_L_sample_freqs_kHz[m] = m * model.getWo() * 4.0 / Math.PI;
        }

        /* clip between peak and peak -50dB, to reduce dynamic range */
        for (m = 1; m <= model.getL(); m++) {
            if (AmdB[m] < (AmdB_peak - 50.0)) {
                AmdB[m] = AmdB_peak - 50.0;
            }
        }

        interp_para(rate_K_vec, 0, rate_L_sample_freqs_kHz, 1, AmdB, 1, model.getL(),
                rate_K_sample_freqs_kHz, 0, NEWAMP1_K);
    }

    private double rate_K_mbest_encode(int[] indexes, double[] x, double[] xq) {
        float[] codebook1 = codebook_newamp.getCodebook(0).getCodeBookArray();
        float[] codebook2 = codebook_newamp.getCodebook(1).getCodeBookArray();
        double[] target = new double[NEWAMP1_K];
        int[] index = new int[MBEST_STAGES];
        double mse, tmp;
        int i, j, n1, n2;

        MBest[] mbest_stage1 = new MBest[MBEST_ENTRIES];
        MBest[] mbest_stage2 = new MBest[MBEST_ENTRIES];

        for (i = 0; i < MBEST_ENTRIES; i++) {
            mbest_stage1[i] = new MBest();
            mbest_stage2[i] = new MBest();
        }

        /* Stage 1 */
        mbest_search(codebook1, x, codebook_newamp.getCodebook(0).getNumberOfElements(),
                mbest_stage1, index);

        /* Stage 2 */
        for (j = 0; j < MBEST_ENTRIES; j++) {
            index[1] = mbest_stage1[j].getIndex(0);

            for (i = 0; i < NEWAMP1_K; i++) {
                target[i] = x[i] - codebook1[NEWAMP1_K * index[1] + i];
            }

            mbest_search(codebook2, target, codebook_newamp.getCodebook(1).getNumberOfElements(),
                    mbest_stage2, index);
        }

        n1 = mbest_stage2[0].getIndex(1);
        n2 = mbest_stage2[0].getIndex(0);

        mse = 0.0;
        for (i = 0; i < NEWAMP1_K; i++) {
            tmp = (double) (codebook1[NEWAMP1_K * n1 + i] + codebook2[NEWAMP1_K * n2 + i]);
            mse += (x[i] - tmp) * (x[i] - tmp);
            xq[i] = tmp;
        }

        indexes[0] = n1;
        indexes[1] = n2;

        return mse;
    }

    private void post_filter_newamp1(double[] vec, double pf_gain) {
        double[] pre = new double[NEWAMP1_K];
        int i;

        /*
         * vec is rate K vector describing spectrum of current frame lets pre-emp
         * before applying PF. 20dB/dec over 300Hz. Postfilter affects energy of
         * frame so we measure energy before and after and normalise. Plenty of room
         * for experiment here as well.
         */
        double e_before = 0.0;
        double e_after = 0.0;

        for (i = 0; i < NEWAMP1_K; i++) {
            pre[i] = 20.0 * Math.log10(rate_K_sample_freqs_kHz[i] / 0.3);
            vec[i] += pre[i];
            e_before += Math.pow(10.0, 2.0 * vec[i] / 20.0);
            vec[i] *= pf_gain;
            e_after += Math.pow(10.0, 2.0 * vec[i] / 20.0);
        }

        double gaindB = 10.0 * Math.log10(e_after / e_before);

        for (i = 0; i < NEWAMP1_K; i++) {
            vec[i] -= gaindB;
            vec[i] -= pre[i];
        }
    }

    private void interp_Wo_v(double[] Wo_, int[] L_, boolean[] voicing_, double Wo1,
            double Wo2, boolean voicing1, boolean voicing2) {
        double c;
        int i;

        // interpolation rate
        for (i = 0; i < 4; i++) {
            voicing_[i] = false;
        }

        if (!voicing1 && !voicing2) {
            for (i = 0; i < 4; i++) {
                Wo_[i] = TAU / 100.0;
            }
        }

        if (voicing1 && !voicing2) {
            Wo_[0] = Wo_[1] = Wo1;
            Wo_[2] = Wo_[3] = TAU / 100.0;
            voicing_[0] = voicing_[1] = true;
        }

        if (!voicing1 && voicing2) {
            Wo_[0] = Wo_[1] = TAU / 100.0;
            Wo_[2] = Wo_[3] = Wo2;
            voicing_[2] = voicing_[3] = true;
        }

        if (voicing1 && voicing2) {
            for (i = 0, c = 1.0; i < 4; i++, c -= 1.0 / 4.0) {
                Wo_[i] = Wo1 * c + Wo2 * (1.0 - c);
                voicing_[i] = true;
            }
        }

        for (i = 0; i < 4; i++) {
            L_[i] = (int) Math.floor(Math.PI / Wo_[i]);
        }
    }

    private void resample_rate_L(Model model, double[][] rate_K_vec, int index) {
        double[] rate_K_vec_term = new double[NEWAMP1_K + 2];
        double[] rate_K_sample_freqs_kHz_term = new double[NEWAMP1_K + 2];
        double[] AmdB = new double[MAX_AMP + 1];
        double[] rate_L_sample_freqs_kHz = new double[MAX_AMP + 1];
        int i;

        /* terminate either end of the rate K vecs with 0dB points */
        rate_K_vec_term[0] = rate_K_vec_term[NEWAMP1_K + 1] = 0.0;
        rate_K_sample_freqs_kHz_term[0] = 0.0;
        rate_K_sample_freqs_kHz_term[NEWAMP1_K + 1] = 4.0;

        for (i = 0; i < NEWAMP1_K; i++) {
            rate_K_vec_term[i + 1] = rate_K_vec[index][i];
            rate_K_sample_freqs_kHz_term[i + 1] = rate_K_sample_freqs_kHz[i];
        }

        for (i = 1; i <= model.getL(); i++) {
            rate_L_sample_freqs_kHz[i] = i * model.getWo() * 4.0 / Math.PI;
        }

        interp_para(AmdB, 1, rate_K_sample_freqs_kHz_term, 0, rate_K_vec_term, 0,
                NEWAMP1_K + 2, rate_L_sample_freqs_kHz, 1, model.getL());

        for (i = 1; i <= model.getL(); i++) {
            model.setA(i, Math.pow(10.0, AmdB[i] / 20.0));
        }
    }

    private void determine_phase(Complex[][] H, int index, Model model) {
        double[] Gdbfk = new double[NS];
        double[] sample_freqs_kHz = new double[NS];
        double[] phase = new double[NS];
        double[] AmdB = new double[MAX_AMP + 1];
        double[] rate_L_sample_freqs_kHz = new double[MAX_AMP + 1];
        int i;

        for (i = 1; i <= model.getL(); i++) {
            AmdB[i] = 20.0 * Math.log10(model.getA(i));
            rate_L_sample_freqs_kHz[i] = (double) i * model.getWo() * 4.0 / Math.PI;
        }

        for (i = 0; i < NS; i++) {
            sample_freqs_kHz[i] = (FS / 1000.0) * (double) i / NEWAMP1_PHASE_NFFT;
        }

        interp_para(Gdbfk, 0, rate_L_sample_freqs_kHz, 1, AmdB, 1, model.getL(), sample_freqs_kHz, 0, NS);
        mag_to_phase(phase, Gdbfk);

        for (i = 1; i <= model.getL(); i++) {
            int temp = (int) Math.floor(0.5 + i * model.getWo() * NEWAMP1_PHASE_NFFT / TAU);
            H[index][i] = new Complex(Math.cos(phase[temp]), Math.sin(phase[temp]));
        }
    }

    private void newamp1_interpolate(double[][] interpolated_surface_, double[] left_vec, double[] right_vec) {
        int i, k;
        double c;

        // (linearly) interpolate 25Hz amplitude vectors back to 100Hz
        for (i = 0, c = 1.0; i < 4; i++, c -= 1.0 / 4.0) {
            for (k = 0; k < NEWAMP1_K; k++) {
                interpolated_surface_[i][k] = left_vec[k] * c + right_vec[k] * (1.0 - c);
            }
        }
    }

    private void newamp1_indexes_to_rate_K_vec(double[] rate_K_vec_, int[] indexes) {
        double[] rate_K_vec_no_mean_ = new double[NEWAMP1_K];
        float[] codebook1 = codebook_newamp.getCodebook(0).getCodeBookArray();
        float[] codebook2 = codebook_newamp.getCodebook(1).getCodeBookArray();
        int k;

        for (k = 0; k < NEWAMP1_K; k++) {
            rate_K_vec_no_mean_[k] = (double) (codebook1[NEWAMP1_K * indexes[0] + k]
                    + codebook2[NEWAMP1_K * indexes[1] + k]);
        }

        post_filter_newamp1(rate_K_vec_no_mean_, 1.5);

        double mean_ = codebook_energy.getCodebook(0).getCodeBookArray(indexes[2]);

        for (k = 0; k < NEWAMP1_K; k++) {
            rate_K_vec_[k] = rate_K_vec_no_mean_[k] + mean_;
        }
    }

    private int quantize(float[] cb, double[] mean, int m) {
        double diff, e;

        int besti = 0;
        double beste = 1E32;

        for (int j = 0; j < m; j++) {
            e = 0.0;

            diff = (double) cb[j] - mean[0];
            e += (diff * diff);

            if (e < beste) {
                beste = e;
                besti = j;
            }
        }

        return besti;
    }

    private int encode_log_Wo(double Wo, int bits) {
        int Wo_levels = 1 << bits;

        double norm = (Math.log10(Wo) - WO_MIN) / (WO_MAX - WO_MIN);
        int index = (int) Math.floor(Wo_levels * norm + 0.5);

        if (index < 0) {
            index = 0;
        }

        if (index > (Wo_levels - 1)) {
            index = Wo_levels - 1;
        }

        return index;
    }

    private double decode_log_Wo(int index, int bits) {
        double step = (WO_MAX - WO_MIN) / (1 << bits);

        return Math.pow(10.0, WO_MIN + step * index);
    }
}
