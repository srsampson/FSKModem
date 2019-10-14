/*
 * Copyright (C) 1993-2017 David Rowe
 *
 * All rights reserved
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package codec2;

public final class Codec2 {
    public static final double TAU = (Math.PI * 2.0);
    //
    public static final int FS = 8000;
    public static final int P_MIN = 20;
    public static final int P_MAX = 160;
    //
    public static final int NEWAMP1_K = 20;
    public static final int FFT_SIZE = 512;
    public static final int NEWAMP1_PHASE_NFFT = 128;
    public static final int NS = (NEWAMP1_PHASE_NFFT / 2 + 1);
    public static final int N_SAMP = 80;
    public static final int M_PITCH = 320;
    public static final int MAX_AMP = 80;
    public static final int NW = 280;
    //
    private final Model encodeModel;
    private final Model[] decodeModels;
    //
    private final double[] prev_rate_K_vec_;
    private final double[] Wo_left;
    private final boolean[] voicing_left;
    //
    private final Codebooknewamp1_energy codebook_energy;
    private final FFT fftencode;
    private final FFT fftdecode;
    private final FFT fftphase;
    private final Pack bitEncode;
    private final Pack bitDecode;
    private final Analyze analyze;
    private final Synthesize synthesize;
    private final Newamp1 newamp;

    public Codec2() {
        fftdecode = new FFT(FFT_SIZE);     // for full duplex
        fftencode = new FFT(FFT_SIZE);     // we need separate FFT's
        fftphase = new FFT(NEWAMP1_PHASE_NFFT);

        bitEncode = new Pack();    // same for pack
        bitDecode = new Pack();    // independant for duplex

        synthesize = new Synthesize(fftdecode);
        analyze = new Analyze(fftencode);
        newamp = new Newamp1(fftphase);
        codebook_energy = new Codebooknewamp1_energy();

        Wo_left = new double[1];
        voicing_left = new boolean[1];
        prev_rate_K_vec_ = new double[NEWAMP1_K];

        /*
         * Pre-allocate the four voice models for decode
         */
        decodeModels = new Model[4];

        for (int i = 0; i < 4; i++) {
            decodeModels[i] = new Model();
        }

        encodeModel = new Model();
    }

    public int codec2_getBitsPerFrame() {
        return 28;
    }

    public int codec2_getBytesPerFrame() {
        return (codec2_getBitsPerFrame() + 7) / 8;
    }

    public int codec2_getSamplesPerFrame() {
        return 320;
    }

    /*
     * The speech output as packed bytes from the 16-bit PCM buffer array
     *
     * The byte array is (codec2_getBitsPerFrame() + 7) / 8 bytes long
     */
    public void codec2_encode(byte[] bits, int index, short[] speech) {
        encode(bits, index, speech);
    }

    /*
     * The speech output as packed bytes from the 16-bit PCM buffer array
     * No offset version
     */
    public void codec2_encode(byte[] bits, short[] speech) {
        encode(bits, 0, speech);
    }

    /*
     * The output of a PCM 16-bit buffer array from an input of packed bytes
     */
    public void codec2_decode(short[] speech, int index, byte[] bits) {
        decode(speech, index, bits);
    }

    /*
     * The output of a PCM 16-bit buffer array from an input of packed bytes
     * No offset version.
     */
    public void codec2_decode(short[] speech, byte[] bits) {
        decode(speech, 0, bits);
    }
    
    /*
     * This provides the user a way to reset everything back to
     * defaults for the current mode.
     */
    public void codec2_reset() {
        bitEncode.reset();
        bitDecode.reset();
        analyze.reset();
        synthesize.reset();
        encodeModel.reset();

        for (int i = 0; i < 4; i++) {
            decodeModels[i].reset();
        }
    }

    /*
     * Encodes 320 speech samples (40ms of speech) into 28 bits. Bits are shifted
     * left into the bytes.
     *
     * In the Coherent Modem, this is called with a 0, or a 320 sample
     * offset index. As the modem processes two codec frames for each modem
     * frame.
     */
    private void encode(byte[] bits, int offset, short[] speech) {
        int[] indexes = new int[4];

        for (int i = 0; i < codec2_getBytesPerFrame(); i++) {
            bits[i] = 0;
        }

        bitEncode.reset();
        
        analyze.analyze_one_frame(encodeModel, speech, offset + 0); 
        analyze.analyze_one_frame(encodeModel, speech, offset + N_SAMP);
        analyze.analyze_one_frame(encodeModel, speech, offset + N_SAMP * 2);
        analyze.analyze_one_frame(encodeModel, speech, offset + N_SAMP * 3);

        newamp.newamp1_model_to_indexes(encodeModel, indexes);

        bitEncode.pack(bits, indexes[0], 9);
        bitEncode.pack(bits, indexes[1], 9);
        bitEncode.pack(bits, indexes[2], 4);
        bitEncode.pack(bits, indexes[3], 6);
    }

    private void decode(short[] speech, int offset, byte[] bits) {
        Complex[][] HH = new Complex[4][MAX_AMP + 1];
        double[][] interpolated_surface_ = new double[4][NEWAMP1_K];
        int[] indexes = new int[4];
        int i, j;

        bitDecode.reset();
        
        indexes[0] = bitDecode.unpack(bits, 9);
        indexes[1] = bitDecode.unpack(bits, 9);
        indexes[2] = bitDecode.unpack(bits, 4);
        indexes[3] = bitDecode.unpack(bits, 6);
        
        for (i = 0; i < 4; i++) {
            for (j = 0; j <= MAX_AMP; j++) {
                HH[i][j] = new Complex();
            }
        }

        newamp.newamp1_indexes_to_model(decodeModels, interpolated_surface_,
                Wo_left, voicing_left, prev_rate_K_vec_, HH, indexes);

        for (i = 0; i < 4; i++) {
            synthesize.synthesize_one_frame(decodeModels[i], speech, offset + N_SAMP * i, HH[i]);
        }
    }

    public double codec2_energy(byte[] bits) {
        int[] indexes = new int[4];

        bitDecode.reset();
        
        indexes[0] = bitDecode.unpack(bits, 9);
        indexes[1] = bitDecode.unpack(bits, 9);
        indexes[2] = bitDecode.unpack(bits, 4);
        indexes[3] = bitDecode.unpack(bits, 6);

        double mean = codebook_energy.getCodebook(0).getCodeBookArray(indexes[2]) - 10.0;

        if (indexes[3] == 0)
            mean -= 10.0;

        return Math.pow(10.0, mean / 10.0);
    }
}
