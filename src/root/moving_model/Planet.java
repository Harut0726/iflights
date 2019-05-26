package root.moving_model;

import root.Vector;
import root.Params;

import java.lang.reflect.Array;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class Planet extends Mover{
    private double orbitalPeriod;
    private int stepOfPartition;

    public Planet(String name, double mass, Vector loc, Vector vel, double orbitalPeriod) {
        this.name = name;
        this.mass = mass;
        this.loc = loc;
        this.vel = vel;
        this.orbitalPeriod = orbitalPeriod;
        this.stepOfPartition = (int) Math.floor(orbitalPeriod / Params.Phys_Params.divParam);
    }

    public Planet(Planet other) {
        this.name = other.name;
        this.mass = other.mass;
        this.loc = new Vector(other.loc.x, other.loc.y);
        this.vel = new Vector(other.vel.x, other.vel.y);
        this.orbitalPeriod = other.orbitalPeriod;
        this.stepOfPartition = other.stepOfPartition;
    }

    public double getOrbitalPeriod() {
        return orbitalPeriod;
    }

    public void setOrbitalPeriod(float orbitalPeriod) {
        this.orbitalPeriod = orbitalPeriod;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getStepOfPartition() {
        return stepOfPartition;
    }

    public void setStepOfPartition(int stepOfPartition) {
        this.stepOfPartition = stepOfPartition;
    }

    public String toString() {
        return String.format("name = %s, mass = %s", name, String.valueOf(mass));
    }

    public static Planet getPlanetByName(ArrayList<Planet> planetsArray, String nameOfPlanet) {
        for(Planet planet: planetsArray){
            if(planet.getName().equals(nameOfPlanet)) {
                return planet;
            }
        }
        throw new RuntimeException("Planet with such name doesn't exits");
    }

    public static void updateVectorsOfPlanets(ArrayList<Planet> planets, ArrayList<Mover> movers) {
        for(Planet planet: planets) {
            if(!planet.getName().equals("Sun"))
                planet.RK4(movers);
        }
    }
}
