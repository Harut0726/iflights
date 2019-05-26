package root.mcts;

import root.Params;

public class State {
    private int daysPassed;
    private double totalDeltaV;

    public int getDaysPassed() {
        return daysPassed;
    }

    public void setDaysPassed(int daysPassed) {
        this.daysPassed = daysPassed;
    }

    public double getTotalDeltaV() {
        return totalDeltaV;
    }

    public void setTotalDeltaV(double totalDeltaV) {
        this.totalDeltaV = totalDeltaV;
    }


    public State(int daysPassed, double totalDeltaV) {
        this.daysPassed = daysPassed;
        this.totalDeltaV = totalDeltaV;
    }

    public State(State other) {
        this.daysPassed = other.daysPassed;
        this.totalDeltaV = other.totalDeltaV;
    }

    public String toString() {
        StringBuilder res = new StringBuilder();
        res.append("daysPassed: " + daysPassed + ", ");
        res.append("totV: " + totalDeltaV);
        return res.toString();
    }

    public boolean isTerminal() {
        int maxDuration = Params.Mission_Params.MISSION_DURATION;
        if(daysPassed > maxDuration) {
            return true;
        }

        double maxV = Params.Mission_Params.MAX_DELTA_V;
        if(totalDeltaV > maxV) {
            return true;
        }

        return false;
    }

    public void updateTotalDeltaV(double deltaV) {
        totalDeltaV += deltaV;
    }
}
