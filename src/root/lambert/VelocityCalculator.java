package root.lambert;

import root.Vector;

public class VelocityCalculator {
    private double v_inf_in_sq, v_inf_out_sq, delta;

    public VelocityCalculator(Vector v_in_inf, Vector v_out_inf) {
        delta = getAngle(v_in_inf, v_out_inf);
        v_inf_in_sq = v_in_inf.magSq();
        v_inf_out_sq = v_out_inf.magSq();
    }

    private double getAngle(Vector v1, Vector v2) {
        return Math.acos((v1.dot(v2)) / (v1.mag() * v2.mag()));
    }

    public double calcDeltaV() {
        double e2_start=0;
        double eps = 0.0001;
        double e2 = newtonRaphson(e2_start, eps);

        double temp = (2 * v_inf_out_sq) / (1-e2);
        double res = Math.sqrt(v_inf_out_sq + temp) - Math.sqrt(v_inf_in_sq + temp);
        return Math.abs(res);
    }

    private double newtonRaphson(double x, double eps) {
        int max_iterations = 10;
        double h = f(x) / dfde2(x);
        int iterations = 0;
        while (Math.abs(h) >= eps)
        {
            h = f(x) / dfde2(x);
            x = x - h;
            iterations++;
            if(iterations > max_iterations) {
                throw new RuntimeException("Too many iterations for Newton's method");
            }
        }
        return x;
    }

    private double f(double e2) {
        return (v_inf_in_sq / v_inf_out_sq * (e2-1) + 1) *
                Math.sin(delta - Math.asin(1/e2)) - 1;
    }

    private double dfde2(double e2) {
        double temp = v_inf_in_sq / v_inf_out_sq;
        return (temp * Math.sin(delta-Math.asin(1/e2))) + (temp * (e2-1) + 1)
                * Math.cos(delta - Math.asin(1/e2)) * 1/(e2*Math.sqrt(e2*e2 - 1));
    }
}
