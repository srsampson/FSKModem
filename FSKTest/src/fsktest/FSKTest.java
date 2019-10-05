/*
 * Copyright (C) 1993-2017 David Rowe, Brady O'Brien
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package fsktest;

import framing.Framing;
import math.Complex;
import fsk.Defines;
import fsk.Demodulator;
import fsk.Modulator;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Random;

public class FSKTest implements Defines {

    private static Random rand;
    private static FileOutputStream dataOut;
    private static FileInputStream dataIn;
    private static Demodulator demod;
    private static Modulator mod;
    private static Framing framer;

    public static void main(String[] args) {
        int lsb, msb, word;

        // write a test file for writing
        try {
            dataOut = new FileOutputStream("/tmp/fskdata.raw");
        } catch (FileNotFoundException | SecurityException e0) {
            System.err.println("Fatal: Couldn't open raw output file: /tmp/fskdata.raw");
            System.exit(-1);
        }

        // Instantiate our modem
        demod = new Demodulator();
        mod = new Modulator();
        framer = new Framing();
        rand = new Random(System.currentTimeMillis());

        // Check to see that the framer works

        // The codec2 data bytes are 28 bits each. That means 4 LSB bits are 0

        byte[] codecData = new byte[]{
            (byte)0xAA, (byte)0xAA, (byte)0x55, (byte)0xF0,
            (byte)0x55, (byte)0x55, (byte)0xAA, (byte)0xF0
        };
        
        byte[] testData = new byte[8];      // This is the codec data decoded from the frame
        byte[] codecBits = new byte[FRAME_BITS];    // The 64-bit Frame Sync + Codec packed

        System.out.print("Here we have 8 bytes of Codec2 Digital Data: ");
        
        for (int i = 0; i < 8; i++) {
            System.out.printf("%02X ", codecData[i]);
        }

        System.out.println("\n\nThis is then packed into 64-bits following the Sync byte:\n");

        framer.getModemFrame(codecBits, codecData);

        for (int i = 0; i < FRAME_BITS; i += 8) {
            for (int j = 0; j < 8; j++) {
                System.out.print(codecBits[i+j] == 1 ? "1" : "0");
            }

            System.out.print(" ");
        }

        framer.getCodecFrame(testData, codecBits);

        System.out.printf("\n\nThis is the verified 8 Codec2 Digital Data bytes: ");

        for (int i = 0; i < 8; i++) {
            System.out.printf("%02X ", testData[i]);
        }

        double snr = demod.fsk_get_snr();
        
        System.out.printf("\n\nFrame SNR: %.02f%n%n", snr);

        double[] modulation = new double[CYCLES * FRAME_SYMBOLS];   // 20 * 32
        byte[] bitpair = new byte[FRAME_BITS];

        try {
            for (int loop = 0; loop < 4; loop++) {  // loop 4 times to give a good quantity of data

                // Modulate 64 bits and write to a file
                for (int i = 0; i < FRAME_BITS; i += 2) {

                    // Select a random FSK tone (1 of 4) or 2-bits
                    int temp = rand.nextInt(4);
                    bitpair[i] = (byte)((temp >>> 1) & 0x1); // MSB
                    bitpair[i + 1] = (byte)(temp & 0x01);    // LSB
                }
                
                mod.modulate(modulation, bitpair);

                // Store the modem cycle audio in little-endian raw signed 16-bit PCM
                for (int l = 0; l < (CYCLES * FRAME_SYMBOLS); l++) {
                    // 50% modulation

                    word = (int) (modulation[l] * 16384.0) & 0xFFFF;
                    lsb = word & 0xFF;
                    msb = (word >>> 8) & 0xFF;
                    dataOut.write(lsb);
                    dataOut.write(msb);
                }

                // Show the original 32 symbols (64 bits) random bits modulated
                for (int i = 0; i < FRAME_BITS; i += 2) {
                    System.out.printf("%d%d", bitpair[i], bitpair[i + 1]);
                }

                System.out.println();
            }
        } catch (IOException e1) {
            System.err.println("Fatal: Couldn't write to audio file");
            System.exit(-1);
        }

        try {
            dataOut.flush();
            dataOut.close();
        } catch (IOException e2) {
        }

        // open the test file for decoding
        try {
            dataIn = new FileInputStream("/tmp/fskdata.raw");
        } catch (FileNotFoundException | SecurityException e3) {
            System.err.println("Fatal: Couldn't open raw audio file: /tmp/fskdata.raw");
            System.exit(-1);
        }

        System.out.println();

        byte[] bits = new byte[FRAME_BITS];              // two bits per symbol (64 bits)
        double[] freqs = new double[NUMBER_OF_TONES];
        Complex[] rx = new Complex[CYCLE_MEMORY];       // 640 + 40

        // Initialize receiver buffer
        for (int i = 0; i < CYCLE_MEMORY; i++) {
            rx[i] = new Complex();
        }

        try {
            for (int loop = 0; loop < 4; loop++) {

                // Now read in the audio for 640 samples and decode the 64-bits
                for (int i = 0; i < FRAME_SAMPLES; i++) {
                    lsb = dataIn.read();
                    msb = dataIn.read();

                    rx[i] = new Complex((double) ((short) (msb << 8 | lsb)) / 16384.0, 0.0);  // +/- 1.0 range
                }

                demod.demodulate(freqs, bits, rx);

                for (int i = 0; i < NUMBER_OF_TONES; i++) {
                    System.out.printf("%.0f Hz%n", freqs[i] * (double) SAMPLERATE); // convert to Hz
                }

                // We have 64-bits now
                for (int i = 0; i < FRAME_BITS; i++) {
                    System.out.print(bits[i]);
                }

                System.out.printf("%nSNR = %.2f%n%n", demod.fsk_get_snr());

            }
        } catch (IOException e4) {
            System.out.println("Decode exception: " + e4.toString());
        }

        try {
            dataIn.close();
        } catch (IOException e5) {
        }
    }
}
