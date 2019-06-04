package root;

public class Params {
    //Mission parameters
    public static class Mission_Params {
        public static final String START_DATE = "1977-08-20"; //15:33 {hh:mm}
        public static final String NAME_OF_LAST_PLANET = "Neptune";
        public static final int MAX_MISSION_DURATION = 4387; //1977-08-20 – 1989-08-24"
        public static final double MAX_DELTA_V = 0.005732894322353824;
    }

    //Physical parameters
    public static class Phys_Params {
        public static final double G_SCALED = 1.4881376098413224e-34; //G * 2.22972472 * pow(10, -24)
        public static final double MU = 0.0002959117; // mass of Sun * G_SCALED
        public static final double DT = 1; //step of integration | days
        public static final double eps = 1e-6;
        public static final int divParam = 29;
    }
}
