package root.mcts;

import java.util.ArrayList;

import root.Params;
import root.Vector;
import root.db.DatabaseUtils;
import root.lambert.LambertSolver;
import root.lambert.VelocityCalculator;
import root.moving_model.Planet;

public class TimeNode extends Node {
    private int flightDuration;
    private Vector v1;
    private Vector v2;

    private Vector planetVelocity;

    public TimeNode(int flightDuration, PlanetNode parent, State state, int level) {
        this.flightDuration = flightDuration;
        this.parent = parent;
        this.state = state;
        this.level = level;
        computationPerformed = false;
    }

    public int getFlightDuration() {
        return flightDuration;
    }

    public void setFlightDuration(int flightDuration) {
        this.flightDuration = flightDuration;
    }

    public Vector getV1() {
        return v1;
    }

    public void setV1(Vector v1) {
        this.v1 = v1;
    }

    public Vector getV2() {
        return v2;
    }

    public void setV2(Vector v2) {
        this.v2 = v2;
    }

    public Vector getPlanetVelocity() {
        return planetVelocity;
    }

    public void setPlanetVelocity(Vector planetVelocity) {
        this.planetVelocity = planetVelocity;
    }

    public String toString() {
        String res = "[duration: " + flightDuration +
                ", reward: " +
                reward +
                "], [state: " +
                state +
                "] " +
                "v1=" +
                v1 +
                "v2=" +
                v2;
        return res;
    }

    @Override
    ArrayList<Node> initChildren(ArrayList<Planet> planetsArr) {
        ArrayList<Node> res = new ArrayList<>();
        for(Planet pl: planetsArr) {
            PlanetNode par = (PlanetNode) parent;
            if(!pl.getName().equals("Sun") && !pl.getName().equals(par.getNameOfPlanet())) {
                PlanetNode planetNode = new PlanetNode(
                        pl.getName(),
                        this,
                        new State(
                                getState().getDaysPassed(),
                                state.getTotalDeltaV()
                        ),
                        level + 1
                );
                res.add(planetNode);
            }
        }
        return res;
    }

    @Override
    public void doComputation() {
        if(!lambertLegComputed()) {
            calcLambertLeg();
        }
        if (!computationPerformed) {
            double deltaV = calcDeltaV();
            state.updateTotalDeltaV(deltaV);
            computationPerformed = true;
        }
    }

    private boolean lambertLegComputed() {
        return v1 != null && v2 != null;
    }

    private void calcLambertLeg() {
        assert !lambertLegComputed():
                "Lambert's leg was already computed";

        PlanetNode planetNodeFrom = getPlanetNodeFrom();
        PlanetNode planetNodeTo = getPlanetNodeTo();

        Vector[] vectorsFrom = DatabaseUtils.getVectorsByNameAndTime(
                "Ephemeris",
                planetNodeFrom.getNameOfPlanet(),
                parent.parent.getState().getDaysPassed()
        );
        Vector[] vectorsTo = DatabaseUtils.getVectorsByNameAndTime(
                "Ephemeris",
                planetNodeTo.getNameOfPlanet(),
                state.getDaysPassed());

        assert state.getDaysPassed() - parent.parent.getState().getDaysPassed() == flightDuration;
        //solving lambert's problem
        Vector[] vels  = LambertSolver.lambertcpp(
                vectorsFrom[0], //r1
                vectorsTo[0],   //r2
                flightDuration,
                Params.Phys_Params.MU
        );

        v1 = vels[0];
        v2 = vels[1];
        planetVelocity = vectorsTo[1];
    }

    private PlanetNode getPlanetNodeFrom() {
        Node dad = parent.parent;
        if(dad.parent == null) {
            return (PlanetNode) dad;
        }
        return (PlanetNode) dad.parent;
    }

    private PlanetNode getPlanetNodeTo() {
        return (PlanetNode) parent;
    }

    private double calcDeltaV() {
        assert lambertLegComputed() &&
                (level==3 || ((TimeNode) parent.parent).lambertLegComputed()):
                "Can't calculate delta v. Lambert's leg isn't computed";
        double res;
        if(level == 3) {
            Vector earthVel = DatabaseUtils.getVectorsByNameAndTime(
                    "Ephemeris",
                    "Earth",
                    0)[1];
            res = Vector.sub(v1, earthVel).mag();
        }
        else {
            Vector plVel = ((TimeNode) parent.parent).getPlanetVelocity();
            Vector vIn = ((TimeNode) parent.parent).getV2();
            Vector vOut = getV1();
            Vector vInRelToPl = Vector.sub(vIn, plVel);
            Vector vOutRelToPl = Vector.sub(vOut, plVel);

            res = new VelocityCalculator(vInRelToPl, vOutRelToPl).calcDeltaV();
        }
        return res;
    }

    @Override
    public boolean isTerminalNode() {
        PlanetNode par = (PlanetNode) parent;
        String nameOfLastPlanet = Params.Mission_Params.NAME_OF_LAST_PLANET;
        return par.getNameOfPlanet().equals(nameOfLastPlanet);
    }
}
