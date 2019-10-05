/*
 * Copyright (C) 1993-2017 David Rowe
 *
 * All rights reserved
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package codec2;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Scanner;

/*
 * I'm using floats here, to save memory space (4 versus 8 bytes each)
 */
public final class Codebooknewamp1 {

    private Scanner scanner;
    private InputStreamReader stream;
    private final Codebook[] vq_cb;
    private final float[] codes0;
    private final float[] codes1;

    /*
     * We have to read these VQ values in, as the number of values is too
     * great for the Java compiler to handle as static constants.
     * We all miss Fortran IV, I'm sure...
     *
     * These text files are set as resources in the JAR file, so the user
     * doesn't have to maintain them. They are in the package 'codebook'
     */
    public Codebooknewamp1() {
        vq_cb = new Codebook[2];

        try {
            stream = new InputStreamReader(getClass().getResourceAsStream("/codebook/train_120_1.txt"));
            scanner = new Scanner(new BufferedReader(stream));
        } catch (Exception fnfe1) {
            System.err.println("Vector file not found '/codebook/train_120_1.txt' " + fnfe1.getMessage());
            System.exit(1);
        }

        int k = scanner.nextInt();
        int m = scanner.nextInt();
        int log2m = (int) (Math.log((double) m) / Math.log(2.0));
        
        codes0 = new float[k * m];
        vq_cb[0] = new Codebook(k, log2m, m, codes0);

        int index = 0;
        while (scanner.hasNextFloat()) {
            codes0[index++] = scanner.nextFloat();
        }

        scanner.close();

        try {
            stream = new InputStreamReader(getClass().getResourceAsStream("/codebook/train_120_2.txt"));
            scanner = new Scanner(new BufferedReader(stream));
        } catch (Exception fnfe2) {
            System.err.println("Vector file not found '/codebook/train_120_2.txt' " + fnfe2.getMessage());
            System.exit(2);
        }

        k = scanner.nextInt();
        m = scanner.nextInt();
        log2m = (int) (Math.log((double) m) / Math.log(2.0));

        codes1 = new float[k * m];
        vq_cb[1] = new Codebook(k, log2m, m, codes1);

        index = 0;
        while (scanner.hasNextFloat()) {
            codes1[index++] = scanner.nextFloat();
        }

        scanner.close();
    }

    public Codebook getCodebook(int index) {
        return vq_cb[index];
    }
}
