/*
 * Copyright (C) 2010 Perens LLC
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package codec2;

public final class Pack {

    private int bitOffset;

    public Pack() {
        bitOffset = 0;
    }

    public void reset() {
        bitOffset = 0;
    }

    protected void pack(byte[] bitArray, int value, int valueBits) {
        do {
            int bI = bitOffset;
            int bitsLeft = 8 - (bI & 0x7);
            int sliceWidth = bitsLeft < valueBits ? bitsLeft : valueBits;
            int wordIndex = bI >>> 3;

            bitArray[wordIndex] |= ((byte) ((value >> (valueBits - sliceWidth)) << (bitsLeft - sliceWidth)));

            bitOffset = bI + sliceWidth;
            valueBits -= sliceWidth;
        } while (valueBits != 0);
    }

    protected int unpack(byte[] bitArray, int valueBits) {
        int field = 0;

        do {
            int bI = bitOffset;
            int bitsLeft = 8 - (bI & 0x7);
            int sliceWidth = bitsLeft < valueBits ? bitsLeft : valueBits;

            field |= (((bitArray[bI >>> 3] >> (bitsLeft - sliceWidth)) & ((1 << sliceWidth) - 1)) << (valueBits - sliceWidth));

            bitOffset = bI + sliceWidth;
            valueBits -= sliceWidth;
        } while (valueBits != 0);

        return field;
    }
}
