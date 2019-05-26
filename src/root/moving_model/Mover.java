package root.moving_model;

import root.Vector;
import java.util.ArrayList;
import root.Params;

public abstract class Mover {
    protected String name;
    protected Vector loc, vel; //location and velocity
    protected double mass;

    public Vector getLoc() {
        return loc;
    }

    public void setLoc(Vector loc) {
        this.loc = loc;
    }

    public Vector getVel() {
        return vel;
    }

    public void setVel(Vector vel) {
        this.vel = vel;
    }

    public double getMass() {
        return mass;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }

    public void RK4(ArrayList<Mover> movers) {
        Vector kv1, kv2, kv3, kv4;
        Vector kr1, kr2, kr3, kr4;

        double h = Params.Phys_Params.DT;

        kr1 = new Vector(vel.x, vel.y);
        kv1 = a(loc, movers);
        kr2 = Vector.add(vel, Vector.mult(kv1, 0.5*h));
        kv2 = a(Vector.add(loc, Vector.mult(kr1, 0.5*h)), movers);
        kr3 = Vector.add(vel, Vector.mult(kv2, 0.5*h));
        kv3 = a(Vector.add(loc, Vector.mult(kr2, 0.5*h)), movers);
        kr4 = Vector.add(vel, Vector.mult(kv3, h));
        kv4 = a(Vector.add(loc, Vector.mult(kr3, h)), movers);

        Vector loc_temp = Vector.add(kr1, Vector.mult(kr2,2));
        loc_temp.add(Vector.mult(kr3,2));
        loc_temp.add(kr4);
        loc_temp.mult(h/6.0);
        loc.add(loc_temp);

        Vector vel_temp = Vector.add(kv1, Vector.mult(kv2,2));
        vel_temp.add(Vector.mult(kv3,2));
        vel_temp.add(kv4);
        vel_temp.mult(h/6.0);
        vel.add(vel_temp);
    }

    public Vector a(Vector loc, ArrayList<Mover> movers) {
        Vector acceleration = new Vector(0,0,0);
        Vector force;

        for(Mover mover: movers) {
            if(!mover.equals(this)) {
                force = attract(mover, loc);
                acceleration.add(force);
            }
        }

        return acceleration;
    }

    public Vector attract(Mover mover, Vector loc) {
        Vector force = Vector.sub(mover.loc, loc);
        double dSquare = force.magSq();
        force.normalize();
        double strength = (Params.Phys_Params.G_SCALED * mover.mass) /
                (dSquare + Params.Phys_Params.eps);
        force.mult(strength);

        return force;
    }

    // this method can be called instead of RK4 for faster (inaccurate) calculations
    public void takeAvgValues(ArrayList<Mover> movers) {
        double h = Params.Phys_Params.DT;

        Vector vel_old = vel;
        Vector a_old = a(loc, movers);

        Vector loc_expect = Vector.add(loc, Vector.mult(vel_old, h));
        Vector a_expect = a(loc_expect, movers);

        Vector vel_temp = Vector.add(a_old, a_expect);
        vel_temp.mult(0.5*h);
        vel = Vector.add(vel_old, vel_temp);

        Vector loc_temp = Vector.add(vel_old, vel);
        loc_temp.mult(0.5*h);
        loc.add(loc_temp);
    }
}
