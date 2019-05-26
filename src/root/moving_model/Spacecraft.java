package root.moving_model;

import root.Vector;

import java.util.ArrayList;

public class Spacecraft extends Mover{
    private double angle;

    public Spacecraft() {}

    public Spacecraft(String name, double mass, Vector loc, Vector vel) {
        this.name = name;
        this.mass = mass;
        this.loc = loc;
        this.vel = vel;
        angle = Math.toRadians(Math.random() * 360);
    }

    public Spacecraft(Spacecraft other) {
        this.name = other.name;
        this.mass = other.mass;
        this.loc = new Vector(other.loc.x, other.loc.y);
        this.vel = new Vector(other.vel.x, other.vel.y);
    }

    public double getAngle() {
        return angle;
    }

    public void setAngle(double angle) {
        this.angle = angle;
    }

    public String toString() {
        return String.format("name = %s, mass = %s", name, String.valueOf(mass));
    }

    public void updateVectors(ArrayList<Mover> movers) {
        RK4(movers);
    }
}
