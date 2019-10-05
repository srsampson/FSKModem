/*
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package math;

public final class Complex {

    private final double real;
    private final double imag;

    public Complex() {
        real = 0.0;
        imag = 0.0;
    }

    public Complex(double re, double im) {
        real = re;
        imag = im;
    }

    public double getReal() {
        return real;
    }

    public double getImaginary() {
        return imag;
    }
}
