package root.lambert;

import root.Vector;
import java.util.Arrays;


public class LambertSolver {
    public static double acosh(double x) {
        return Math.log(x + Math.sqrt(x*x-1));
    }

    public static double asinh(double x) {
        return Math.log(x + Math.sqrt(x*x+1));
    }

    public static double sinh(double x) {
        return 0.5 * (Math.exp(x)-Math.exp(-x));
    }

    public static int sign(double x){
        if(x > 0) {
            return 1;
        }
        if(x < 0) {
            return -1;
        }
        return 0;
    }

    static double x2tof(double x, double s, double c, int lw) {
        double am, a, alfa, beta;

        am = s / 2;
        a = am / (1 - x * x);
        if (x < 1)
        {
            beta = 2 * Math.asin(Math.sqrt((s - c) / (2 * a)));
            if (lw != 0) beta = -beta;
            alfa = 2 * Math.acos(x);
        }
        else
        {
            alfa = 2 * acosh(x);
            beta = 2 * asinh(Math.sqrt((s - c) / (-2 * a)));
            if (lw != 0) beta = -beta;
        }

        if (a > 0)
            return (a * Math.sqrt(a)* ((alfa - Math.sin(alfa)) - (beta - Math.sin(beta))));
        else
            return (-a * Math.sqrt(-a)*((sinh(alfa) - alfa) - (sinh(beta) - beta)));
    }

    static void vett(double[] vet1, double[] vet2, double[] prod) {
        prod[0] = (vet1[1] * vet2[2] - vet1[2] * vet2[1]);
        prod[1] = (vet1[2] * vet2[0] - vet1[0] * vet2[2]);
        prod[2] = (vet1[0] * vet2[1] - vet1[1] * vet2[0]);
    }

    static void vers(double[] V_in, double[] Ver_out) {
        double v_mod = 0;
        int i;

        for (i = 0; i < 3; i++)
        {
            v_mod += V_in[i] * V_in[i];
        }

        double sqrtv_mod = Math.sqrt(v_mod);

        for (i = 0; i < 3; i++)
        {
            Ver_out[i] = V_in[i] / sqrtv_mod;
        }
    }

    //Izzo's method c++
    public static Vector[] lambertcpp(Vector r1vec, Vector r2vec, double t, double mu) {
        double[] v1 = new double[3], v2 = new double[3];
        double a, p, theta;
        int iter;
        int lw = 0;

        double V, T,
                r2_mod = 0,
                dot_prod = 0,
                c,
                s,
                am,
                lambda,
                x, x1, x2, y1, y2, x_new = 0, y_new, err,
                alfa, beta, psi, eta, eta2, sigma1, vr1,
                vt1, vt2, vr2, R=0.0;
        int i_count, i;
        double tolerance = 1e-11;
        double[] r1 = new double[3],
                r2 = new double[3],
                r2_vers = new double[3],
                ih_dum = new double[3],
                ih = new double[3],
                dum = new double[3];

        assert t > 0 : "t must be > 0";

        double[] r1_in = {r1vec.x, r1vec.y, 0};
        double[] r2_in = {r2vec.x, r2vec.y, 0};
        for (i = 0; i < 3; i++) {
            r1[i] = r1_in[i];
            r2[i] = r2_in[i];
            R += r1[i] * r1[i];
        }

        R = Math.sqrt(R);
        V = Math.sqrt(mu / R);
        T = R / V;

        t /= T;
        for (i = 0; i < 3; i++) {
            r1[i] /= R;
            r2[i] /= R;
            r2_mod += r2[i] * r2[i];
        }

        r2_mod = Math.sqrt(r2_mod);

        for (i = 0; i < 3; i++)
            dot_prod += (r1[i] * r2[i]);

        theta = Math.acos(dot_prod / r2_mod);

        c = Math.sqrt(1 + r2_mod * (r2_mod - 2.0 * Math.cos(theta)));
        s = (1 + r2_mod + c) / 2.0;
        am = s / 2.0;
        lambda = Math.sqrt(r2_mod) * Math.cos(theta / 2.0) / s;
        x1 = Math.log(0.4767);
        x2 = Math.log(1.5233);
        y1 = Math.log(x2tof(-.5233, s, c, lw)) - Math.log(t);
        y2 = Math.log(x2tof(.5233, s, c, lw)) - Math.log(t);

        err = 1;
        i_count = 0;
        while ((err>tolerance) && (y1 != y2)) {
            i_count++;
            x_new = (x1*y2 - y1*x2) / (y2 - y1);
            y_new = Math.log(x2tof(Math.exp(x_new) - 1, s, c, lw)) - Math.log(t);
            x1 = x2;
            y1 = y2;
            x2 = x_new;
            y2 = y_new;
            err = Math.abs(x1 - x_new);
        }
        iter = i_count;
        x = Math.exp(x_new) - 1;

        a = am / (1 - x*x);

        if (x < 1) {
            beta = 2 * Math.asin(Math.sqrt((s - c) / (2*a)));
            if (lw != 0) beta = -beta;
            alfa = 2 * Math.acos(x);
            psi = (alfa - beta) / 2;
            eta2 = 2 * a*Math.pow(Math.sin(psi), 2) / s;
            eta = Math.sqrt(eta2);
        }
        else {
            beta = 2 * asinh(Math.sqrt((c - s) / (2 * a)));
            if (lw != 0) beta = -beta;
            alfa = 2 * acosh(x);
            psi = (alfa - beta) / 2;
            eta2 = -2 * a * Math.pow(sinh(psi), 2) / s;
            eta = Math.sqrt(eta2);
        }

        p = (r2_mod / (am * eta2)) * Math.pow(Math.sin(theta / 2), 2);
        sigma1 = (1 / (eta * Math.sqrt(am)))* (2 * lambda * am - (lambda + x * eta));
        vett(r1, r2, ih_dum);
        vers(ih_dum, ih);

        if (lw != 0) {
            for (i = 0; i < 2; i++)
                ih[i] = -ih[i];
        }

        vr1 = sigma1;
        vt1 = Math.sqrt(p);
        vett(ih, r1, dum);

        for (i = 0; i < 3; i++)
            v1[i] = vr1 * r1[i] + vt1 * dum[i];

        vt2 = vt1 / r2_mod;
        vr2 = -vr1 + (vt1 - vt2) / Math.tan(theta / 2);
        vers(r2, r2_vers);
        vett(ih, r2_vers, dum);
        for (i = 0; i < 3; i++)
            v2[i] = vr2 * r2[i] / r2_mod + vt2 * dum[i];

        for (i = 0; i < 3; i++)
        {
            v1[i] *= V;
            v2[i] *= V;
        }
        a *= R;
        p *= R;
        Vector[] res = {
                new Vector(v1[0], v1[1]),
                new Vector(v2[0], v2[1])
        };
        return res;
    }


    //Izzo's method
    static Vector[] lambert1(Vector r1vec, Vector r2vec, double tf, double m, double muC) {
        double tol = 1e-14;
        boolean bad = false;
        int days = 86400;

        double r1 = r1vec.mag();
        r1vec.div(r1);
        double V = Math.sqrt(muC/r1);
        r2vec.div(r1);
        double T = r1/V;
        tf = tf*days/T;

        double mr2vec = r2vec.mag();
        Matrix r1matr = new Matrix(new double[]{r1vec.x, r1vec.y});
        Matrix r2matr = new Matrix(new double[]{r2vec.x, r2vec.y});
        Matrix r1r2 = r1matr.multiply(r2matr.transpose());
        double r1r2Num = r1r2.toDouble();
        double dth = Math.acos(Math.max(-1, Math.min(1, (r1r2Num)/mr2vec)));

        int leftbranch = sign(m);
        int longway = sign(tf);
        m = Math.abs(m);
        tf = Math.abs(tf);
        if (longway < 0)  dth = 2*Math.PI - dth;

        double c = Math.sqrt(1 + Math.pow(mr2vec, 2) - 2*mr2vec*Math.cos(dth));
        double s = (1 + mr2vec + c)/2;
        double a_min = s/2;
        double Lambda = Math.sqrt(mr2vec)*Math.cos(dth/2)/2;
        Vector crossprd = r1vec.cross(r2vec);
        double mcr = crossprd.mag();
        Vector nrmunit = Vector.div(crossprd, (float)mcr);

        double logt = Math.log(tf);

        double inn1, inn2, x1, x2;
        if (m == 0) {
            inn1 = -0.5233;
            inn2 = 0.5233;
            x1 = Math.log(1 + inn1);
            x2 = Math.log(1 + inn2);
        }
        else {
            if (leftbranch < 0) {
                inn1 = -0.5234;
                inn2 = 0.5234;
            }
            else {
                inn1 = 0.7234;
                inn2 = 0.5234;
            }
            x1 = Math.tan(inn1*Math.PI/2);
            x2 = Math.tan(inn2*Math.PI/2);
        }

        double[] xx = {inn1, inn2};
        double[] aa = {
                a_min / (1 - inn1*inn1),
                a_min / (1 - inn2*inn2)
        };
        double[] bbeta = {
                longway * 2*Math.asin(Math.sqrt((s-c)/2/aa[0])),
                longway * 2*Math.asin(Math.sqrt((s-c)/2/aa[1]))
        };
        //not sure START min
        double[] tempArr = {1, inn1, inn2};
        Arrays.sort(tempArr);
        double aalfa = 2*Math.acos(Math.max(-1, tempArr[0]));
        //not sure END min

        double[] y12 = {
                aa[0]*Math.sqrt(aa[0])*
                        ((aalfa - Math.sin(aalfa)) -
                                (bbeta[0]-Math.sin(bbeta[0])) + 2*Math.PI*m),
                aa[1]*Math.sqrt(aa[1])*
                        ((aalfa - Math.sin(aalfa)) -
                                (bbeta[1]-Math.sin(bbeta[1])) + 2*Math.PI*m),
        };

        double y1, y2;
        if (m == 0) {
            y1 = Math.log(y12[0]) - logt;
            y2 = Math.log(y12[1]) - logt;
        }
        else {
            y1 = y12[0] - tf;
            y2 = y12[1] - tf;
        }

        double err = 1000000;
        int iterations = 0;
        double xnew = 0;
        while(err > tol) {
            iterations += 1;
            xnew = (x1*y2 - y1*x2) / (y2-y1);
            double x;
            if (m == 0)
                x = Math.exp(xnew) - 1;
            else
                x = Math.atan(xnew)*2/Math.PI;
            double a = a_min / (1-x*x);
            double beta, alfa;
            if (x < 1) {
                beta = longway * 2*Math.asin(Math.sqrt((s-c)/2/a));
                alfa = 2*Math.acos(Math.max(-1, Math.min(1, x)));
            }
            else {
                alfa = (float)2*acosh(x);
                beta = longway * 2*asinh(Math.sqrt((s-c)/(-2*a)));
            }
            double tof;
            if (a > 0)
                tof = a*Math.sqrt(a)*((alfa - Math.sin(alfa)) -
                        (beta-Math.sin(beta)) + 2*Math.PI*m);
            else
                tof = -a*Math.sqrt(-a)*((sinh(alfa) - alfa) -
                        (sinh(beta) - beta));

            double ynew;
            if (m == 0)
                ynew = Math.log(tof);
            else
                ynew = tof - tf;

            x1 = x2;
            x2 = xnew;
            y1 = y2;
            y2 = ynew;
            err = Math.abs(x1 - xnew);

            if(iterations > 15)
                bad = true;
            break;
        }

        if (bad == true)
            return lambert2(
                    Vector.mult(r1vec, r1),
                    Vector.mult(r2vec, r1),
                    longway*tf*T,
                    leftbranch*m,
                    muC);

        double x;
        if (m == 0)
            x = Math.exp(xnew) - 1;
        else
            x = Math.atan(xnew)*2/Math.PI;

        double a = a_min/(1-x*x);

        double beta, alfa, psi, eta2, eta;
        if (x < 1) {
            beta = longway * 2*Math.asin(Math.sqrt((s-c)/2/a));
            alfa = 2*Math.acos(Math.max(-1, Math.min(1, x)));
            psi = (alfa-beta)/2;
            //not sure START sin^2
            eta2 = 2*a*Math.pow(Math.sin(psi), 2)/s;
            //not sure END sin^2
            eta = Math.sqrt(eta2);
        }
        else {
            beta = longway * 2*asinh(Math.sqrt((c-s)/2/a));
            alfa = 2*acosh(x);
            psi = (alfa-beta)/2;
            //not sure START sinh^2
            eta2 = -2*a*Math.pow(sinh(psi), 2)/s;
            //not sure END sinh^2
            eta = Math.sqrt(eta2);
        }

        Vector ih = Vector.mult(nrmunit, longway);

        Vector r2n = Vector.div(r2vec, mr2vec);

        Vector crsprd1 = ih.cross(r1vec);
        Vector crsprd2 = ih.cross(r2n);

        double Vr1 = 1/eta/Math.sqrt(a_min) * (2*Lambda*a_min - Lambda - x*eta);
        //not sure START sin^2
        double Vt1 = Math.sqrt(mr2vec/a_min/eta2 * Math.pow(Math.sin(dth/2), 2));
        //not sure END sin^2
        double Vt2 = Vt1/mr2vec;
        double Vr2 = (Vt1 - Vt2)/Math.tan(dth/2) - Vr1;

        Vector V1 = Vector.add(
                Vector.mult(r1vec, Vr1),
                Vector.mult(crsprd1, Vt1)
        );
        V1.mult(V);
        Vector V2 = Vector.add(
                Vector.mult(r2n, (float)Vr2),
                Vector.mult(crsprd2, (float)Vt2)
        );
        V2.mult((float)V);

        return new Vector[] {V1, V2};
    }


    //Gooding's method
    static Vector[] lambert2(Vector r1vec, Vector r2vec, double tf, double m, double muC) {
        double x0;

        double tol = 1e-12;

        double r1 = r1vec.mag();
        double r2 = r2vec.mag();

        Vector r1unit = Vector.div(r1vec, (float)r1);
        Vector r2unit = Vector.div(r2vec, (float)r2);
        Vector crsprod = r1vec.cross(r2vec);
        double mcrsprd = crsprod.mag();
        Vector th1unit = Vector.div(crsprod, mcrsprd).cross(r1unit);
        Vector th2unit = Vector.div(crsprod, mcrsprd).cross(r2unit);

        Matrix r1matr = new Matrix(new double[]{r1vec.x, r1vec.y});
        Matrix r2matr = new Matrix(new double[]{r2vec.x, r2vec.y});
        Matrix r1r2 = r1matr.multiply(r2matr.transpose());
        double r1r2Num = r1r2.toDouble();
        double dth = Math.acos(Math.max(-1, Math.min(1, r1r2Num/r1/r2)));

        int longway = sign(tf);
        tf = Math.abs(tf);
        if(longway < 0) {
            dth = dth - 2*Math.PI;
        }

        int leftbranch = sign(m);
        m = Math.abs(m);

        double c = Math.sqrt(r1*r1 + r2*r2 - 2*r1*r2*Math.cos(dth));
        double s = (r1 + r2 + c) / 2;
        double T = Math.sqrt(8*muC/Math.pow(s, 3)) * tf;
        double q = Math.sqrt(r1*r2)/s * Math.cos(dth/2);

        //not sure START Lancaster output
        double[] lancParse = LancasterBlanchard(0, q, m);
        double T0 = lancParse[0];
        //not sure END Lancaster output
        double Td = T0 - T;

        //not sure START mod
        double phr = (2*Math.atan2(1 - Math.pow(q, 2), 2*q)) % (2*Math.PI);
        //not sure END mod

        Matrix V1 = new Matrix(1, 3);
        Matrix V2 = V1;
        //extremal_distances = [NaN, NaN]

        if (m == 0) {
            double x01 = T0*Td/4/T;
            double x02, x03;
            if (Td > 0)
                x0 = x01;
            else {
                x01 = Td/(4 - Td);
                x02 = -Math.sqrt(-Td/(T+T0/2));
                double W = x01 + 1.7*Math.sqrt(2 - phr/Math.PI);
                if(W >= 0)
                    x03 = x01;
                else
                    //not sure START .^.*
                    x03 = x01 + Math.pow(-W, 1/16)*(x02 - x01);
                //not sure END .^.*
                double lambda = 1 + x03*(1 + x01)/2 -
                        0.03*Math.pow(x03, 2)*Math.sqrt(1 + x01);
                x0 = lambda*x03;
            }
            if(x0 < -1)
                throw new RuntimeException("No solution");
        }
        else {
            double xM0;
            double xMpi = 4/(3*Math.PI*(2*m + 1));
            if (phr < Math.PI)
                xM0 = xMpi*Math.pow(phr/Math.PI, 1/8);
            else if (phr > Math.PI)
                xM0 = xMpi*(2 - Math.pow(2-phr/Math.PI, 1/8));
            else
                xM0 = 0;

            double xM = xM0;
            //not sure START inf
            double Tp = 100000;
            //not sure END inf
            int iterations = 0;

            while(Math.abs(Tp) > tol) {
                iterations += 1;
                double[] derivatives = LancasterBlanchard(xM, q, m);
                Tp = derivatives[1];
                double Tpp = derivatives[2];
                double Tppp = derivatives[3];

                double xMp = xM;
                xM = xM - 2*Tp*Tpp / (2*Math.pow(Tpp, 2) - Tp*Tppp);
                //not sure START mod
                if(iterations % 7 != 0)
                    xM = (xMp+xM)/2;
                //not sure END mod
                if(iterations > 25)
                    throw new RuntimeException("exitflag = -2");
            }

            if(xM < -1 || xM > 1)
                throw new RuntimeException("exitflag = -1");

            double[] correspondingTime;
            correspondingTime = LancasterBlanchard(xM, q, m);
            double TM = correspondingTime[0];

            if (TM > T)
                throw new RuntimeException("T should lie above the minimum T");

            double TmTM = T - TM;
            double T0mTM = T0 - TM;
            double[] forParsing = LancasterBlanchard(xM, q, m);
            Tp = forParsing[1];
            double Tpp = forParsing[2];


            double W;
            double x03;
            double x00;
            if(leftbranch > 0) {
                double x = Math.sqrt(TmTM / (Tpp/2 + TmTM/Math.pow(1-xM, 2)));
                W = xM + x;
                W = 4*W/(4 + TmTM) + Math.pow(1-W, 2);
                x0 = x*(1 - (1 + m + (dth - 1/2)) /
                        (1 + 0.15*m)*x*(W/2 + 0.03*x*Math.sqrt(W))) + xM;

                if (x0 > 1)
                    throw new RuntimeException("exitflag = -1");
            }
            else {
                if(Td > 0)
                    x0 = xM - Math.sqrt(TM/(Tpp/2 - TmTM*(Tpp/2/T0mTM - 1/Math.pow(xM, 2))));
                else {
                    x00 = Td / (4 - Td);
                    W = x00 +1.7*Math.sqrt(2*(1 - phr));
                    if (W >= 0)
                        x03 = x00;
                    else
                        x03 = x00 - Math.sqrt(Math.pow(-W, 1/8))*(x00+Math.sqrt(-Td/(1.5*T0 - Td)));
                    W = 4/(4 - Td);
                    double lambda = (1 +(1 + m + 0.24*(dth - 1/2)) /
                            (1 + 0.15*m)*x03*(W/2 - 0.03*x03*Math.sqrt(W)));
                    x0 = x03*lambda;
                }

                if (x0 < -1)
                    throw new RuntimeException("estimate might not give solutions");
            }
        }

        double x = x0;
        double Tx = 1000000;
        int iterations = 0;

        while(Math.abs(Tx) > tol) {
            iterations += 1;

            double[] toParse = LancasterBlanchard(x, q, m);
            Tx = toParse[0];
            double Tp = toParse[1];
            double Tpp = toParse[2];

            Tx = Tx - T;
            double xp = x;
            x = x - 2*Tx*Tp / (2*Math.pow(Tp, 2) - Tx*Tpp);
            if (iterations % 7 != 0)
                x = (xp+x)/2;

            if(iterations > 25)
                throw new RuntimeException("exitflag = -2");
        }

        double gamma = Math.sqrt(muC*s/2);
        double sigma, rho;
        double z;
        if (c == 0) {
            sigma = 1;
            rho = 0;
            z = Math.abs(x);
        }
        else {
            sigma = 2*Math.sqrt(r1*r2/Math.pow(c, 2)) * Math.sin(dth/2);
            rho = (r1 - r2) / c;
            z = Math.sqrt(1 + Math.pow(q, 2)*(Math.pow(x, 2) -1));
        }

        double Vr1 = gamma*((q*z - x) - rho*(q*z + x)) / r1;
        Vector Vr1vec = Vector.mult(r1unit, Vr1);
        double Vr2 = -gamma*((q*z - x) + rho*(q*z + x)) / r2;
        Vector Vr2vec = Vector.mult(r2unit, Vr2);

        double Vtan1 = sigma * gamma * (z + q*x) / r1;
        Vector Vtan1vec = Vector.mult(th1unit, (float)Vtan1);
        double Vtan2 = sigma * gamma * (z + q*x) / r2;
        Vector Vtan2vec = Vector.mult(th2unit, (float)Vtan2);

        Vector V1ans = Vector.add(Vtan1vec, Vr1vec);
        Vector V2ans = Vector.add(Vtan2vec, Vr2vec);
        return new Vector[] {V1ans, V2ans};
    }


    static double[] LancasterBlanchard(double x, double q, double m) { ;
        if(x < -1)
            x = Math.abs(x) - 2;
        else if (x == -1)
            x = x + 2.2204e-16;

        double E = x*x - 1;

        double T, Tp, Tpp, Tppp;

        if(x == 1) {
            T = 4/3*(1-Math.pow(q, 3));
            Tp = 4/5 * (Math.pow(q, 5) - 1);
            Tpp = Tp + 120/70 * (1 - Math.pow(q, 7));
            Tppp = 3*(Tpp - Tp) + 2400/1080 * (Math.pow(q, 9) - 1);
        }
        else if(Math.abs(x-1) < 1e-2) {
            double[] sigParse1 = sigmax(-E);
            double sig1 = sigParse1[0];
            double dsigdx1 = sigParse1[1];
            double d2sigdx21 = sigParse1[2];
            double d3sigdx31 = sigParse1[3];

            double[] sigParse2 = sigmax(-E*q*q);
            double sig2 = sigParse2[0];
            double dsigdx2 = sigParse2[1];
            double d2sigdx22 = sigParse2[2];
            double d3sigdx32 = sigParse2[3];

            T = sig1 - Math.pow(q, 3) * sig2;
            Tp = 2*x*(Math.pow(q, 5)*dsigdx2 - dsigdx1);
            Tpp = Tp/x + 4*Math.pow(x,2)*(d2sigdx21 - Math.pow(q, 7)*d2sigdx22);
            Tppp = 3*(Tpp-Tp/x)/x + 8*x*x*(Math.pow(q, 9)*d3sigdx32 - d3sigdx31);
        }
        else {
            double y = Math.sqrt(Math.abs(E));
            double z = Math.sqrt(1 + Math.pow(q, 2)*E);
            double f = y*(z - q*x);
            double g = x*z - q*E;

            //float factor1 = E < 0 ? 1: 0;
            //float factor2 = E > 0 ? 1: 0;
            //float d = factor1 * (atan2(f, g) + PI*m) +
            //  factor2 * log(max(0, f+g));
            double d;
            if (E<0)
                d = Math.atan2(f, g) + Math.PI*m;
            else if (E == 0)
                d = 0;
            else
                d = Math.log(Math.max(0, f+g));

            T = 2*(x - q*z - d/y)/E;
            Tp = (4 - 4*Math.pow(q, 3)*x/z - 3*x*T)/E;
            Tpp = (-4*Math.pow(q,3)/z * (1-Math.pow(q,2)*Math.pow(x,2)/Math.pow(z,2)) -
                    3*T - 3*x*Tp)/E;
            Tppp = (4*Math.pow(q, 3)/Math.pow(z,2) * ((1-Math.pow(q,2)*Math.pow(x,2)/Math.pow(z,2)) +
                    2*Math.pow(q,2)*x/Math.pow(z,2)*(z-x)) - 8*Tp - 7*x*Tpp)/ E;
        }
        double[] res = {T, Tp, Tpp, Tppp};
        return res;
    }

    static double[] sigmax(double y) {
        double[] anArr = new double[] {
                4.000000000000000e-001,
                2.142857142857143e-001,
                4.629629629629630e-002,
                6.628787878787879e-003,
                7.211538461538461e-004,
                6.365740740740740e-005,
                4.741479925303455e-006,
                3.059406328320802e-007,
                1.742836409255060e-008,
                8.892477331109578e-010,
                4.110111531986532e-011,
                1.736709384841458e-012,
                6.759767240041426e-014,
                2.439123386614026e-015,
                8.203411614538007e-017,
                2.583771576869575e-018,
                7.652331327976716e-020,
                2.138860629743989e-021,
                5.659959451165552e-023,
                1.422104833817366e-024,
                3.401398483272306e-026,
                7.762544304774155e-028,
                1.693916882090479e-029,
                3.541295006766860e-031,
                7.105336187804402e-033
        };

        Matrix an = new Matrix(anArr);
        an = an.transpose();
        assert an.M == 25 && an.N == 1;

        double[] powersArr = new double[25];
        for(int i = 0; i < 25; i++) {
            powersArr[i] = Math.pow(y, i+1);
        }

        Matrix powers = new Matrix(powersArr); // y; y^2;..;y^25


        Matrix sig = powers.multiply(an).plus(4/3);
        //Matrix sig = new Matrix(1,3);

        //not sure START dsigma/dx
        double[] temp = new double[25]; //1 2y .. 25y^24
        temp[0] = 1;
        for (int i=1; i < 25; i++)
            temp[i] = (i+1) * powersArr[i-1];

        Matrix dsigdx = new Matrix(temp).multiply(an);
        //not sure END dsigma/dx

        //not sure START d2sigma/dx2
        double[] temp2 = new double[25]; //0 2 6y .. 25*24y^23
        temp2[0] = 0;
        temp2[1] = 2;
        for (int i=2; i < 25; i++)
            temp2[i] = (i+1)*i * powersArr[i-2];

        Matrix d2sigdx2 = new Matrix(temp2).multiply(an);
        //not sure END d2sigma/dx2

        //not sure START d3sigma/dx3
        double[] temp3 = new double[25]; //0 0 6 24y .. 25*24*23y^22
        temp3[0] = 0;
        temp3[1] = 0;
        temp3[2] = 6;
        for (int i=3; i < 25; i++)
            temp3[i] = (i+1)*i*(i-1) * powersArr[i-3];

        Matrix d3sigdx3 = new Matrix(temp3).multiply(an);
        //not sure END d3sigma/dx3

        double[] res = {sig.toDouble(), dsigdx.toDouble(),
                d2sigdx2.toDouble(), d3sigdx3.toDouble()};
        return res;
    }
}
