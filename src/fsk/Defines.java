/*
 * Copyright (C) 1993-2017 David Rowe, Brady O'Brien
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package fsk;

public interface Defines {

    double TAU = (Math.PI * 2.0);
    double RATE = 25.0;                                       // Vocoder Rate
    double PERIOD = (1.0 / RATE);                             // 40 ms
    //
    int FFTSIZE = 512;
    int SYMBOLRATE = 400;
    int FIRST_TONE = 800;
    int TONE_SEPARATION = SYMBOLRATE;
    int NUMBER_OF_TONES = 4;
    int SAMPLERATE = 8000;
    int CYCLES = SAMPLERATE / SYMBOLRATE;                     // 20
    int SYMBOLS = (int) ((double) SYMBOLRATE * PERIOD);       // 16 Symbols (32 bits)
    int OVERSAMPLE = 10;
    int NSTASH = CYCLES * 4;                                  // 80
    int EST_SPACE = (SYMBOLRATE - (SYMBOLRATE / 5));          // 320
    int EST_MIN = (SYMBOLRATE / 2);
    int EST_MAX = (SAMPLERATE / 2) - SYMBOLRATE;
    int FRAME_SYMBOLS = (SYMBOLS * 2);
    int FRAME_SAMPLES = (FRAME_SYMBOLS * CYCLES);
    int FRAME_BITS = (FRAME_SYMBOLS * 2);                     // 64
    int CYCLE_MEMORY = (FRAME_SAMPLES + (CYCLES * 2));        // 640 + 40
    //
    double TC = 0.95 * (double) FFTSIZE / (double) SAMPLERATE;// .0608
}
