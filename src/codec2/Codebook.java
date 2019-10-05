/*
 * Copyright (C) 1993-2017 David Rowe
 *
 * All rights reserved
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package codec2;

public final class Codebook {

    private final int k;         // dimension of the vector
    private final int log2m;     // log base 2 of size m
    private final int m;         // number of vector elements
    private final float[] cb;    // the actual codebook array

    protected Codebook(int kval, int log2mval, int mval, float[] cbval) {
        k = kval;
        log2m = log2mval;
        m = mval;
        cb = cbval;
    }

    protected int getDimension() {
        return k;
    }

    protected int getLogBase2Elements() {
        return log2m;
    }

    protected int getNumberOfElements() {
        return m;
    }

    protected float[] getCodeBookArray() {
        return cb;
    }

    protected float getCodeBookArray(int index) {
        return cb[index];
    }
}
