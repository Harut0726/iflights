package root.lambert;

import root.Params;
import root.TimeFormatter;
import root.Vector;

import java.sql.Time;

public class VoyagerImpulseCalculator {
    private static Encounter
            earth = new Encounter(
                        new Vector(8.561576983560533E-01, -5.389692490908495E-01),
                        new Vector(8.881530388789968E-03, 1.450102105675779E-02),
                 "1977-08-20"
            ),
            jupiter = new Encounter(
                        new Vector(-3.936629795576152E+00, 3.593219485219530E+00),
                        new Vector(-5.181517737944351E-03, -5.229916901145773E-03),
                    "1979-07-09"
            ),
            saturn = new Encounter(
                        new Vector(-9.392457004359585E+00, -1.951806812557789E+00),
                        new Vector(8.393586571245021E-04, -5.480879611705111E-03),
                    "1981-08-25"
            ),
            uranus = new Encounter(
                        new Vector(-3.667338021225589E+00, -1.876748128811125E+01),
                        new Vector(3.835797772464359E-03, -9.329136564981506E-04),
                    "1986-01-24"
            ),
            neptune = new Encounter(
                        new Vector(5.999621390193559E+00, -2.960814343253063E+01),
                        new Vector(3.051416457593553E-03, 6.422718862963416E-04),
                    "1989-08-24"
            );

    public static double VoyagerImpulse() {
        Vector[] vels1 = getVels(earth, jupiter);
        double deltaV1 = Vector.sub(vels1[0], earth.v).mag();

        Vector[] vels2 = getVels(jupiter, saturn);
        double deltaV2 = (new VelocityCalculator(vels1[1], vels2[0])).calcDeltaV();

        Vector[] vels3 = getVels(saturn, uranus);
        double deltaV3 = (new VelocityCalculator(vels2[1], vels3[0])).calcDeltaV();

        Vector[] vels4 = getVels(uranus, neptune);
        double deltaV4 = (new VelocityCalculator(vels3[1], vels4[0])).calcDeltaV();

//        Vector[] vels5 = getVels(earth, neptune);
//        double deltaV5 = Vector.sub(vels5[0], earth.v).mag();

        return deltaV1 + deltaV2 + deltaV3 + deltaV4;
    }

    private static Vector[] getVels(Encounter encounter1, Encounter encounter2) {
        double mu = Params.Phys_Params.MU;
        long t = TimeFormatter.daysBetweenDates(encounter1.date, encounter2.date);
        return LambertSolver.lambertcpp(encounter1.r, encounter2.r, t, mu);
    }

    public static double convertToKmPerSec(double v) {
        double au = 149_597_871; //km
        double day = 86400; //seconds
        return v * au / day;
    }
}

