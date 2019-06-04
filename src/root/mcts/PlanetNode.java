package root.mcts;

import java.util.ArrayList;

import root.Params;
import root.moving_model.Planet;


public class PlanetNode extends Node {
    private String nameOfPlanet;

    public PlanetNode(String nameOfPlanet, Node parent, State state, int level) {
        this.nameOfPlanet = nameOfPlanet;
        this.parent = parent;
        this.state = state;
        this.level = level;
    }

    public PlanetNode(String nameOfPlanet){
        this.nameOfPlanet = nameOfPlanet;
    }

    public String getNameOfPlanet() {
        return nameOfPlanet;
    }

    public void setNameOfPlanet(String nameOfPlanet) {
        this.nameOfPlanet = nameOfPlanet;
    }

    public String toString() {
        String res = "[name: " + nameOfPlanet + "], " +
                "[visits: " +
                nVisits +
                ", reward: " +
                reward + "], [state: " +
                state +
                "]";
        return res;
    }

    //@Override
    ArrayList<Node> initChildren2(ArrayList<Planet> planetsArr) {
        ArrayList<Node> res = new ArrayList<>();
        int step = getStepByName(this.nameOfPlanet, planetsArr);
        ArrayList<Integer> times = divideTime(state.getDaysPassed(), step);
        for(int time: times) {
            TimeNode timeNode = new TimeNode(
                    time - state.getDaysPassed(),
                    this,
                    new State(time+state.getDaysPassed(), state.getTotalDeltaV()),
                    level+1
            );
            res.add(timeNode);
        }
        return res;
    }

    @Override
    ArrayList<Node> initChildren(ArrayList<Planet> planetsArr) {
        ArrayList<Node> res = new ArrayList<>();
        int step = getStepByName(nameOfPlanet, planetsArr);
        int days_passed = state.getDaysPassed();
        int duration = Params.Mission_Params.MAX_MISSION_DURATION;
        int st = days_passed / step;

        for (int k = st+1; k*step <= duration; ++k) {
            int fl_dur = k*step - days_passed;

            res.add(new TimeNode(
                        fl_dur,
                        this,
                        new State(k*step, state.getTotalDeltaV()),
                        level+1
                    )
            );
        }

        return res;
    }

    private int getStepByName(String nameOfPl, ArrayList<Planet> planetsArr) {
        Planet pl = Planet.getPlanetByName(planetsArr, nameOfPl);
        return pl.getStepOfPartition();
    }

    private ArrayList<Integer> divideTime(int daysPassed, int step) {
        ArrayList<Integer> res = new ArrayList<>();
        int duration = Params.Mission_Params.MAX_MISSION_DURATION;

        for(int days=step; days < duration; days += step) {
            if(days > daysPassed) {
                res.add(days);
            }
        }
        return res;
    }

    @Override
    public void doComputation() {

    }

    @Override
    public boolean isTerminalNode() {
        return false;
    }

}
