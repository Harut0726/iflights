package root;

import root.db.DatabaseUtils;
import root.lambert.VoyagerImpulseCalculator;
import root.mcts.PlanetNode;
import root.mcts.State;
import root.mcts.Node;
import root.mcts.Tree;
import root.moving_model.ModelLoader;
import root.moving_model.Planet;

import java.util.ArrayList;


public class MainClass {
    public static void main(String[] args) {
        String path = "Ephemeris_1977_(15-33).txt";
        DatabaseUtils.connectToDB("src/root/db/solarDB.db");

        ArrayList<Planet> planets = ModelLoader.loadPlanetsArray(path);

        Tree tree = new Tree();

        Node cur = tree.getRoot();

        int nSimulations = 50000;
        System.out.println("Number of simulations: " + nSimulations);

        while (!cur.isTerminalNode() && !cur.getState().isTerminal()) {
            for (int i = 0; i < nSimulations; i++) {
                cur.selectAction(planets);
            }
            cur = cur.selectChildWithMaxReward();
            if(cur == null) {
                System.out.println("Insufficient number of simulations");
            }
            System.out.println("SELECTED: " + cur);
        }
    }
}
