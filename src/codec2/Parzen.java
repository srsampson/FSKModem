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

/*
 * This whole class should be garbage collected after run,
 * since it will never be used again
 */
public final class Parzen {
    
    private static final int TW = 40;   // Trapezoidal synthesis window overlap
    
    public Parzen(double[] Pn) {        // Parzen window
        double win;
        int i;
        
        // Generate Parzen window in time domain
        
        for (i = 0; i < N_SAMP / 2 - TW; i++) {
            Pn[i] = 0.0;
        }

        win = 0.0;
        
        for (i = N_SAMP / 2 - TW; i < N_SAMP / 2 + TW; win += 1.0 / (2 * TW), i++) {
            Pn[i] = win;
        }

        for (i = N_SAMP / 2 + TW; i < 3 * N_SAMP / 2 - TW; i++) {
            Pn[i] = 1.0;
        }

        win = 1.0;
        
        for (i = 3 * N_SAMP / 2 - TW; i < 3 * N_SAMP / 2 + TW; win -= 1.0 / (2 * TW), i++) {
            Pn[i] = win;
        }

        for (i = 3 * N_SAMP / 2 + TW; i < N_SAMP * 2; i++) {
            Pn[i] = 0.0;
        }
    }
}
