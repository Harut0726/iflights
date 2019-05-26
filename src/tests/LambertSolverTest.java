package tests;

import root.Params;
import root.TimeFormatter;
import root.Vector;
import root.lambert.LambertSolver;

import static org.junit.Assert.*;

public class LambertSolverTest {

    @org.junit.Test
    public void lambertcpp() {
        double MU = 0.0002959117;
        Vector jupiterLocation = new Vector(-3.93662E+00, 3.59321E+00);
        float tof_Jupiter = (float) TimeFormatter.daysBetweenDates(
                Params.Mission_Params.START_DATE,
                "1979-07-09");
        Vector[] vels = LambertSolver.lambertcpp(
                new Vector(8.562066294346069E-01, -5.389742916708254E-01),
                jupiterLocation,
                tof_Jupiter,
                MU);

        double eps = 0.0001;
        assertTrue(Math.abs(vels[0].x - 0.010250) < eps);
        assertTrue(Math.abs(vels[0].y - 0.019996) < eps);
        assertTrue(Math.abs(vels[1].x - (-0.005520)) < eps);
        assertTrue(Math.abs(vels[1].y - (-0.000714)) < eps);
    }
}

