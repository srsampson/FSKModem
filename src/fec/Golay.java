/*
 * Public Domain (p) 2016, Tomas Härdin
 */
package fec;

/**
 * Golay Forward Error Correction (FEC) [23,12,7]
 *
 * Encoder/Decoder
 *
 * Golay is only used on the data frames, not speech frames
 */
public class Golay {

    // since we want to avoid bit-reversing we bit-reverse the polynomial instead
    private static final int POLY = 0x00000C75;            // 0xAE3 reversed
    //private static final int ALT_POLY =  0x00000AE3;     // 0xC75 reversed
    private static final int X11 = (1 << 11);

    public int golayEncode(int c) {
        return ((syndrome(c << 11) << 12) | c) & 0x7FFFFF;
    }

    public int golayEncodeWithParity(int c) {
        // parity is computed on all bits in the perfect G23 code
        int ret = golayEncode(c);
        int p = ret;

        p ^= (p >>> 16);
        p ^= (p >>> 8);
        p ^= (p >>> 4);
        p ^= (p >>> 2);
        p ^= (p >>> 1);

        return (ret | ((p & 0x1) << 23)) & 0x7FFFFF;
    }

    /**
     * Given a 23 bit codeword, return the 12 bit corrected value
     * 
     * @param c an int representing the 23 bit codeword
     * @param ret an int pointer to the corrected value
     * @return true on error, false on success
     */
    public boolean golayDecode(int c, int[] ret) {
        int s, c2;

        c = unrotate(c, 12);

        for (int x = 0; x < 23; x++) {
            s = syndrome(c);

            if (popcount(s) <= 3) {
                ret[0] = unrotate(c ^ s, x) & 0xFFF;
                return false;
            }

            for (int t = 0; t < 23; t++) {
                c2 = c ^ (1 << t);

                s = syndrome(c2);

                if (popcount(s) <= 2) {
                    ret[0] = unrotate(c2 ^ s, x) & 0xFFF;
                    return false;
                }
            }

            //rotate
            c = (c >>> 1) | ((c & 1) << 22);
        }

        // return true on error
        return true;
    }

    private int syndrome(int c) {
        for (int x = 11; x >= 0; x--) {
            if ((c & (X11 << x)) != 0) {
                c ^= (POLY << x);
            }
        }

        return c & 0x7FFFFF;
    }

    private int popcount(int c) {
        int ret = 0;

        while (c != 0) {
            if ((c & 1) == 1) {
                ret++;
            }

            c >>>= 1;
        }

        return ret;
    }

    /**
     * Unrotate the codeword value
     * 
     * @param c an int representing the codeword
     * @param x an int representing the number of shifts
     * @return an int representing the result
     */
    private int unrotate(int c, int x) {
        int val1 = (c << x);
        int val2 = (val1 >>> (23 - x));

        return (val1 | val2) & 0x7FFFFF;
    }
}
