/*
 * See Knuth TAOCP vol 2, 3rd edition, page 232
 */
package fsk;

public final class Statistics {

    private double oldM;
    private double newM;
    private double oldS;
    private double newS;
    private int n;

    public Statistics() {
        n = 0;
        oldM = 0.0;
        newM = 0.0;
        oldS = 0.0;
        newS = 0.0;
    }

    public void Clear() {
        n = 0;
    }

    public void Push(double x) {
        n++;

        if (n == 1) {
            oldM = newM = x;
            oldS = 0.0;
        } else {
            newM = oldM + (x - oldM) / n;
            newS = oldS + (x - oldM) * (x - newM);

            // set up for next iteration
            oldM = newM;
            oldS = newS;
        }
    }

    public int NumDataValues() {
        return n;
    }

    public double Mean() {
        return (n > 0) ? newM : 0.0;
    }

    public double Variance() {
        return ((n > 1) ? newS / (n - 1) : 0.0);
    }

    public double StandardDeviation() {
        double v = Variance();

        return (v == 0.0) ? 0.0 : Math.sqrt(v);
    }
}
