package root.lambert;

import root.Vector;

public class Encounter {
    public Vector r, v;
    public String date;

    public Encounter(Vector r, Vector v, String date) {
        this.r = r;
        this.v = v;
        this.date = date;
    }
}
