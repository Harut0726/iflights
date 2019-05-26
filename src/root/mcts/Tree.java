package root.mcts;

import java.util.ArrayList;
import java.util.Arrays;

public class Tree {
    private Node root;

    public Tree() {
        root = new PlanetNode(
                "Earth",
                null,
                new State(0, 0),
                1
        );
        root.setnVisits(0);
        root.setReward(0);
        root.setChildren(initChildren());
    }

    private ArrayList<Node> initChildren() {
        ArrayList<Node> res = new ArrayList<>(Arrays.asList(
                new PlanetNode("Mercury"),
                new PlanetNode("Venus"),
                new PlanetNode("Mars"),
                new PlanetNode("Jupiter"),
                new PlanetNode("Saturn"),
                new PlanetNode("Uranus"),
                new PlanetNode("Neptune")
        ));
        for(Node tn: res) {
            tn.setnVisits(0);
            tn.setReward(0);
            tn.setParent(root);
            tn.setState(root.getState());
            tn.setLevel(2);
        }
        return res;
    }

    public Node getRoot() {
        return root;
    }

    public void setRoot(Node root) {
        this.root = root;
    }
}
