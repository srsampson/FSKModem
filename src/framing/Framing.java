/*
 * Copyright (C) 1993-2017 David Rowe, Brady O'Brien
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package framing;

public final class Framing {

    public static final int FRAME_BITS = 64;
    public static final int DATA_BITS = 56;
    public static final int SYNC_BITS = 8;
    public static final int CODEC_BITS = 28;
    public static final int CODEC_BYTES = 8;

    private static final byte[] SYNCWORD = new byte[]{0, 1, 1, 0, 0, 1, 1, 1};  // 0x67
    private static final byte[] DATAWORD = new byte[]{1, 1, 1, 1, 0, 0, 1, 0};  // 0xF2

    /*
     * Data Format Definition
     *
     * Voice1[00:07]
     * Voice1[08:15]
     * Voice1[16:23]
     * Voice1[24:27] + Voice2[00:03]
     * Voice2[04:11]
     * Voice2[12:19]
     * Voice2[20:27]
     */
    private boolean frameSync;          // state is sync or not sync
    //
    private int frameIndex;             // index into circular bit buffer
    private int frameMissedCount;       // How many Sync have been missed
    private int framesSinceLastSync;    // How many bits since the last Sync
    //
    private long frameTotalSyncBits;    // Total RX-ed bits of SyncWord
    private long frameTotalSyncErrors;  // Total errors in UW bits
    //
    private final byte[] frameBits = new byte[FRAME_BITS];

    public Framing() {
        frameSync = false;

        frameIndex = 0;
        framesSinceLastSync = 0;
        frameMissedCount = 0;

        frameTotalSyncBits = 0L;
        frameTotalSyncErrors = 0L;
    }

    // Get a single bit out of an MSB-first packed byte array
    private byte unpackMSB(byte[] data, int offset, int index) {
        return (byte) (data[offset + (index >>> 3)] >>> (7 - (index & 0x7)) & 0x1);
    }

    // See if the syncWord is where it should be, to within a tolerance
    private boolean matchSyncWord(byte[] bits, int tolerance, int[] diffCount) {
        int diff = 0;

        diffCount[0] = 0;

        // Start bit pointer where Sync should be
        int ibit = frameIndex;

        if (ibit >= FRAME_BITS) {
            ibit -= FRAME_BITS;
        }

        // Walk through and match bits in frame with bits of syncWord

        for (int i = 0; i < SYNC_BITS; i++) {
            if (bits[ibit] != SYNCWORD[i]) {
                diff++;
            }

            ibit++;

            if (ibit >= FRAME_BITS) {
                ibit = 0;
            }
        }

        diffCount[0] = diff;

        return (diff <= tolerance);
    }

    private void extractCodecFrame(byte[] codecBytes, byte[] bits) {
        int i;

        // Initialize the 8 bytes
        for (i = 0; i < CODEC_BYTES; i++) {
            codecBytes[i] = 0;
        }

        // Skip past the Sync Word
        int index = frameIndex + SYNC_BITS;

        if (index >= FRAME_BITS) {
            index -= FRAME_BITS;
        }

        // Extract and pack first voice codec 28-bit frame, MSB first

        for (i = 0; i < CODEC_BITS; i++) {

            codecBytes[i >>> 3] |= (bits[index++] & 0x1) << (7 - (i & 0x7));

            if (index >= FRAME_BITS) {
                index = 0;
            }
        }

        index = frameIndex + SYNC_BITS + CODEC_BITS;

        if (index >= FRAME_BITS) {
            index -= FRAME_BITS;
        }

        // Extract and pack second voice codec 28-bit frame, MSB first

        for (i = 0; i < CODEC_BITS; i++) {

            codecBytes[4 + (i >>> 3)] |= (bits[index++] & 0x1) << (7 - (i & 0x7));

            if (index >= FRAME_BITS) {
                index = 0;
            }
        }
    }

    /**
     * Method to return the 64-bit modem data frame
     *
     * This modem data frame contains the Sync Word and two Codec2 700C Frames
     *
     * @param bits is the bit array containing the modem frame bits
     * @param codecBytes is the 8-bytes of codec2 frame bytes
     */
    public void getModemFrame(byte[] bits, byte[] codecBytes) {
        int i, ibit;

        // Fill out frame with prototype
        for (i = 0; i < SYNC_BITS; i++) {
            bits[i] = SYNCWORD[i];
        }

        // Fill out first 28-bit codec2 700C block
        ibit = 0;

        for (i = 0; i < CODEC_BITS; i++) {
            bits[i + SYNC_BITS] = unpackMSB(codecBytes, 0, ibit);
            ibit++;
        }

        // Fill out second 28-bit codec2 700C block
        ibit = 0;

        for (i = 0; i < CODEC_BITS; i++) {
            bits[i + (SYNC_BITS + CODEC_BITS)] = unpackMSB(codecBytes, 4, ibit);
            ibit++;
        }
    }

    /**
     * Method to return the codec2 frame bytes given a 64-bit modem frame
     *
     * Locate the Sync Word if able, and extract the two Codec2 700C Frames
     *
     * @param codecBytes is the 8-bytes of codec2 frame bytes
     * @param modemFrameBits a 64-bit modem frame
     * @return boolean true if frame has been extracted successfully
     */
    public boolean getCodecFrame(byte[] codecBytes, byte[] modemFrameBits) {
        int[] diffCount = new int[1];

        boolean frameExtracted = false;

        for (int i = 0; i < FRAME_BITS; i++) {

            frameBits[frameIndex++] = modemFrameBits[i];

            if (frameIndex >= FRAME_BITS) {
                frameIndex -= FRAME_BITS;
            }

            if (frameSync == true) {
                // Already synchronized, just wait till syncWord is back where it should be

                framesSinceLastSync++;

                // syncWord should be here. We're sunk, so deframe anyway
                if (framesSinceLastSync == FRAME_BITS) {
                    framesSinceLastSync = 0;

                    if (!matchSyncWord(frameBits, 1, diffCount)) {
                        frameMissedCount++;
                    } else {
                        frameMissedCount = 0;
                    }

                    // If we go over the miss tolerance, go into no-sync
                    if (frameMissedCount > 3) {
                        frameSync = false;
                    }

                    extractCodecFrame(codecBytes, frameBits);
                    frameExtracted = true;

                    frameTotalSyncBits += (long) SYNC_BITS;
                    frameTotalSyncErrors += (long) diffCount[0];
                }
            } else if (matchSyncWord(frameBits, 0, diffCount)) {
                // found sync

                frameSync = true;
                framesSinceLastSync = 0;
                frameMissedCount = 0;

                extractCodecFrame(codecBytes, frameBits);
                frameExtracted = true;

                frameTotalSyncBits += (long) SYNC_BITS;
                frameTotalSyncErrors += (long) diffCount[0];
            }
        }

        return frameExtracted;
    }

    public int getFrameSizeInBits() {
        return FRAME_BITS;
    }

    public int getVoiceCodecSizeInBytes() {
        return CODEC_BYTES;   // 4 bytes per codec2 700 bit/s frame
    }

    public boolean IsFrameSync() {
        return frameSync;
    }

    public long getFrameTotalSyncBits() {
        return frameTotalSyncBits;
    }

    public long getFrameTotalSyncErrors() {
        return frameTotalSyncErrors;
    }
}
