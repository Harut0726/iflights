package root.moving_model;

import com.google.gson.Gson;
import com.google.gson.JsonObject;
import root.Vector;

import javax.imageio.IIOException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

public class ModelLoader {
    private static final Gson gson = new Gson();

    private static class Entity {
        String name;
        double mass;
        double rx, ry, vx, vy;
        double orbitalPeriod;


        private Planet convertToPlanet() {
            return new Planet(name, mass , new Vector(rx, ry), new Vector(vx, vy), orbitalPeriod);
        }

        private Spacecraft convertToSpacecraft() {
            return new Spacecraft(name, mass, new Vector(rx, ry), new Vector(vx, vy));
        }
    }

    public static ArrayList<Planet> loadPlanetsArray(String pathToFile) {
        ArrayList<Planet> res = new ArrayList<>();
        String jsonArray = readFile(pathToFile);
        Entity[] entities = gson.fromJson(jsonArray, Entity[].class);
        for(Entity entity: entities) {
            res.add(entity.convertToPlanet());
        }
        return res;
    }

    public static Spacecraft loadSpacecraft(String pathToFile) {
        String json = readFile(pathToFile);
        return gson.fromJson(json, Entity.class).convertToSpacecraft();
    }

    private static String readFile(String pathToFile) {
        String res = null;
        try {
            byte[] encoded = Files.readAllBytes(Paths.get(pathToFile));
            res =  new String(encoded);
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }
}
