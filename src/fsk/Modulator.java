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

public final class Modulator implements Defines {

    private final Complex[] oscillator;
    private Complex phase;

    public Modulator() {
        oscillator = new Complex[NUMBER_OF_TONES];
        phase = ComplexMath.cexp(new Complex());

        // Create all the oscillator tones
        
        for (int i = 0, tone = FIRST_TONE; i < NUMBER_OF_TONES; i++, tone += TONE_SEPARATION) {
            oscillator[i] = ComplexMath.cexp(new Complex(0.0, TAU * ((double) tone / (double) SAMPLERATE)));
        }
    }

    public void modulate(double[] baseband, byte[] bits) {
        int index = 0;

        for (int i = 0; i < FRAME_BITS; i += 2) {
            int tone = (int)((bits[i] << 1 | bits[i + 1])) & 0x3;
            Complex dph = oscillator[tone];

            // Output signal based on the tone selected
            for (int j = 0; j < CYCLES; j++) {
                phase = ComplexMath.times(phase, dph);
                baseband[index + j] = phase.getReal();
            }

            // Normalize TX phase to prevent drift
            phase = ComplexMath.normalize(phase);

            index += CYCLES;
        }
    }
}
