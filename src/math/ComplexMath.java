/*
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package math;

public final class ComplexMath {

    private ComplexMath() {
    }

    public static Complex add(Complex a, Complex b) {
        return new Complex(a.getReal() + b.getReal(), a.getImaginary() + b.getImaginary());
    }

    public static Complex minus(Complex a, Complex b) {
        return new Complex(a.getReal() - b.getReal(), a.getImaginary() - b.getImaginary());
    }

    public static Complex times(Complex a, Complex b) {
        return new Complex(a.getReal() * b.getReal() - a.getImaginary() * b.getImaginary(), a.getReal() * b.getImaginary() + a.getImaginary() * b.getReal());
    }

    public static Complex times(Complex a, double alpha) {
        return new Complex(a.getReal() * alpha, a.getImaginary() * alpha);
    }

    public static Complex timesConjugate(Complex a, Complex b) {
        return new Complex(a.getReal() * b.getReal() + a.getImaginary() * b.getImaginary(), a.getImaginary() * b.getReal() - a.getReal() * b.getImaginary());
    }

    public static Complex conjugate(Complex a) {
        return new Complex(a.getReal(), -a.getImaginary());
    }

    public static Complex divide(Complex a, double b) {
        return new Complex(a.getReal() / b, a.getImaginary() / b);
    }

    public Complex divide(Complex a, Complex b) {
        double m = b.getReal() * b.getReal() + b.getImaginary() * b.getImaginary();
        return new Complex((a.getReal() * b.getReal() + a.getImaginary() * b.getImaginary()) / m, (a.getImaginary() * b.getReal() - a.getReal() * b.getImaginary()) / m);
    }

    public static Complex negate(Complex a) {
        return new Complex(-a.getReal(), -a.getImaginary());
    }

    public static double square(Complex a) {
        return (a.getReal() * a.getReal()) + (a.getImaginary() * a.getImaginary());
    }

    public static double absolute(Complex a) {
        return Math.hypot(a.getReal(), a.getImaginary());
    }

    public static Complex cexp(Complex a) {
        if (a.getReal() == 0.0) {
            return new Complex(Math.cos(a.getImaginary()), Math.sin(a.getImaginary()));
        } else {
            return new Complex(Math.exp(a.getReal()) * Math.cos(a.getImaginary()), Math.exp(a.getReal()) * Math.sin(a.getImaginary()));
        }
    }
    
    public static double carg(Complex a) {
        return Math.atan2(a.getImaginary(), a.getReal());
    }
 
    public static Complex normalize(Complex a) {
        double mag = absolute(a);
        return new Complex(a.getReal() / mag, a.getImaginary() / mag);
    }
}
