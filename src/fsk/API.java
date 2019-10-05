/*
 * Copyright (C) 1993-2017 David Rowe, Brady O'Brien
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package fsk;

import codec2.Codec2;

import framing.Framing;
import math.Complex;

public final class API {

    private final int speech_samples;
    private final int samples_per_codec_frame;

    private double squelch_threshold;
    private boolean squelch_enabled;

    private final Modulator modulator;
    private final Demodulator demodulator;
    private final Framing framing;
    private final Codec2 codec2;

    public API() {
        modulator = new Modulator();
        demodulator = new Demodulator();
        framing = new Framing();
        codec2 = new Codec2();
        
        squelch_threshold = -100.0;
        squelch_enabled = false;

        samples_per_codec_frame = codec2.codec2_getSamplesPerFrame();   // 320
        speech_samples = 2 * samples_per_codec_frame;                   // 640
    }

    public void send(double[] fdm, short[] speech_in) {
        byte[] bits = new byte[Defines.FRAME_BITS];
        byte[] codec_bytes = new byte[4];
        byte[] cframe = new byte[8];
        int i, j;

        // Send the audio frame pair to the Vocoder to convert to digital
        
        for (j = 0; j < 2; j++) {
            codec2.codec2_encode(codec_bytes, j * 320, speech_in);

            for (i = 0; i < 4; i++) {
                cframe[(j * 4) + i] = codec_bytes[i];
            }
        }

        // Frame the codec bytes and add Sync word

        framing.getModemFrame(bits, cframe);
        
        // create the modulator waveform

        modulator.modulate(fdm, bits);
    }

    /**
     * Utility method to access the codec bytes before vocoding
     * 
     * @param codec_bytes the bytes to be sent to the vocoder
     * @param signal the complex signal containing Sync and Data Frames
     */
    public void receiveCodecBytes(byte[] codec_bytes, Complex[] signal) {
        byte[] frame_bits = new byte[Defines.FRAME_BITS];
        double[] freqs = new double[4];

        demodulator.demodulate(freqs, frame_bits, signal);
        framing.getCodecFrame(codec_bytes, frame_bits);
    }

    /**
     * Method to process the received complex data frame
     * 
     * @param speech_out the 640 sample PCM audio decoded from the codec
     * @param signal the complex signal containing Sync and Data Frames
     * @return the number of PCM audio samples (640)
     */
    public int receive(short[] speech_out, Complex[] signal) {
        byte[] frame_bits = new byte[Defines.FRAME_BITS];
        byte[] modem_bytes = new byte[8];
        byte[] codec_bytes = new byte[4];
        double[] freqs = new double[4];
        int i, j, k;

        // Find the Sync bytes, and pull out the codec bytes into bits
        
        demodulator.demodulate(freqs, frame_bits, signal);

        // Find the Sync Word and put the bits into bytes
        // 4 bytes of 28-bits for each of the two frames

        framing.getCodecFrame(modem_bytes, frame_bits);

        if (framing.IsFrameSync() == true) {

            // Each modem frame contains two codec frames
            
            for (j = 0; j < 2; j++) {
                
                // Get each of the four bytes

                for (k = 0; k < 4; k++) {
                    codec_bytes[k] = modem_bytes[j * 4 + k];
                }
                
                // convert the digital 28-bit vocoder data to PCM audio
                
                codec2.codec2_decode(speech_out, j * 320, codec_bytes);

                if ((squelch_enabled == true) && (demodulator.fsk_get_snr() < squelch_threshold)) {
                    // get rid of both pairs
                    
                    for (i = 0; i < speech_samples; i++) {
                        speech_out[i] = (short) 0;
                    }
                    
                    return speech_samples;
                }
            }

            return speech_samples;
        } else {
            // no sync
            
            if (squelch_enabled == true) {
                // squelch enabled
                
                for (i = 0; i < speech_samples; i++) {
                    speech_out[i] = (short) 0;
                }
            } else {
                // squelch disabled
                
                for (i = 0; i < speech_samples; i++) {
                    speech_out[i] = (short) signal[i].getReal();
                }
            }
        
            return speech_samples;
        }
    }

    /**
     * get the FSK Modem status
     * 
     * @param sync a pointer to a boolean representing the modem sync state
     * @param snr a pointer to a double representing the modem SNR in dB
     */
    public void getModemStatistics(boolean[] sync, double[] snr) {
        sync[0] = framing.IsFrameSync();
        snr[0] = demodulator.fsk_get_snr();
    }

    /**
     * Get the number of FSK frame PCM samples available
     * 
     * @return an int representing the number of FSK frame PCM samples available
     */
    public int getNIN() {
        return demodulator.fsk_get_nin();
    }

    /**
     * Set the squelch on or off
     * 
     * @param val a boolean representing the squelch being enabled
     */
    public void setSquelchEnable(boolean val) {
        squelch_enabled = val;
    }

    /**
     * Set the threshold where squelch will activate (dB)
     * 
     * @param val a double representing the threshold in dB for squelch activation
     */
    public void setSquelchThreshold(double val) {
        squelch_threshold = val;
    }

    /**
     * Get the number of PCM speech samples in a codec frame (320)
     * 
     * @return an int representing the number of PCM speech samples in a codec frame
     */
    public int getSamplesperCodecFrame() {
        return samples_per_codec_frame;
    }

    /**
     * Get the number of PCM speech samples in an FSK frame (640)
     * 
     * @return an int representing the number of PCM speech samples in an FSK frame
     */
    public int getSpeechSamples() {
        return speech_samples;
    }

    /**
     * Get the Maximum number of FSK Frame Samples (680)
     * 
     * @return an int representing the maximum size of a frame
     */
    public int getMaximumModemSamples() {
        return Defines.CYCLE_MEMORY;
    }

    /**
     * Get the Nominal number of FSK Frame Samples (640)
     * 
     * @return an int representing the nominal size of a frame
     */
    public int getNominalModemSamples() {
        return Defines.FRAME_SAMPLES;
    }
}