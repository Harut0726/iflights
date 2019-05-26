package root.mcts;

import java.util.*;

import root.Params;
import root.moving_model.Planet;

public abstract class Node {
    //Type type;
    protected State state;
    protected Node parent;
    protected ArrayList<Node> children;
    //int nActions = 5;
    protected Random rnd = new Random();
    protected double epsilon = 1e-6;
    protected double nVisits, reward;
    public boolean computationPerformed = false;
    protected int level;

    abstract ArrayList<Node> initChildren(ArrayList<Planet> planetsArr);
    public abstract void doComputation();
    public abstract boolean isTerminalNode();

    public State getState() {
        return state;
    }

    public void setState(State state) {
        this.state = state;
    }

    public Node getParent() {
        return parent;
    }

    public void setParent(Node parent) {
        this.parent = parent;
    }

    public double getnVisits() {
        return nVisits;
    }

    public void setnVisits(double nVisits) {
        this.nVisits = nVisits;
    }

    public double getReward() {
        return reward;
    }

    public void setReward(double reward) {
        this.reward = reward;
    }

    public ArrayList<Node> getChildren() {
        return children;
    }

    public void setChildren(ArrayList<Node> children) {
        this.children = children;
    }

    public int getLevel(){
        return level;
    }

    public void setLevel(int level) {
        this.level = level;
    }

    public void selectAction(ArrayList<Planet> planetsArr) {
        List<Node> visited = new LinkedList<>();
        Node cur = this;
        cur.doComputation();
        visited.add(this);
        while (!cur.isLeaf()) {
            cur = cur.select();
            cur.doComputation();
            visited.add(cur);
        }
        cur.expandNode(planetsArr);
        Node newNode = cur.select();
        if(newNode == null) {
            backPropagation(visited, 0);
            return;
        }
        visited.add(newNode);
        newNode.doComputation();
        double reward = simulateFlight(newNode, planetsArr);
        backPropagation(visited, reward);
    }

    private boolean isLeaf() {
        if(children == null) {
            return true;
        }
        return children.size() == 0;
    }

    //UCT
    public Node select() {
        Node selected = null;
        double bestValue = Double.MIN_VALUE;
        for (Node child : children) {
            double uctValue = child.reward + Math.sqrt(2*Math.log(nVisits+1) /
                    (child.nVisits + epsilon)) + rnd.nextDouble() * epsilon;
            if (uctValue > bestValue) {
                selected = child;
                bestValue = uctValue;
            }
        }
        return selected;
    }

    public Node selectChildWithMaxReward() {
        double maxReward = Double.MIN_VALUE;
        Node maxNode = null;
        for(Node node: children) {
            if(node.reward > maxReward) {
                maxReward = node.reward;
                maxNode = node;
            }
        }
        return maxNode;
    }

    private void expandNode(ArrayList<Planet> planetsArr) {
        children = initChildren(planetsArr);
        for(Node child: children) {
            child.setnVisits(0);
            child.setReward(0);
        }
    }

    private double simulateFlight(Node nodeFrom, ArrayList<Planet> planetsArr) {
        Node cur = nodeFrom;

        while(!cur.isTerminalNode()) {
            if(cur.state.isTerminal()) {
                return 0;
            }
            cur.doComputation();
            ArrayList<Node> children = cur.initChildren(planetsArr);
            if(children.size() == 0) {
                return 0;
            }
            cur = pickRandomChild(children);
        }
        cur.doComputation();
        double reward = cur.calcReward();
        return reward;
    }

    private Node pickRandomChild(ArrayList<Node> children) {
        assert (children.size() > 0): "Length of children is 0";
        return children.get(rnd.nextInt(children.size()));
    }

    private double calcReward() {
        double totDeltaV = state.getTotalDeltaV();
        double maxDeltaV = Params.Mission_Params.MAX_DELTA_V;
        return Math.max(0, 1 - (totDeltaV / maxDeltaV));
    }

    private void backPropagation(List<Node> visitedArr, double simulationResult) {
        for(Node node: visitedArr) {
            node.updateStats(simulationResult);
        }
    }

    private void updateStats(double reward) {
        nVisits++;
        if(this.reward < reward) {
            this.reward = reward;
        }
    }
}
