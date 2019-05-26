package root.db;

import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Connection;
import java.sql.ResultSet;
import root.moving_model.Planet;
import root.Vector;


public class DatabaseUtils {
    private static Connection conn;
    private static PreparedStatement pstmt;
    private static ResultSet resSet;


    public static void connectToDB(String path) {
        try {
            conn = null;
            String url = "jdbc:sqlite:" + path;
            conn = DriverManager.getConnection(url);
        }
        catch (SQLException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void insert(String tableName, EphemerisData ep) {
        String query = String.format("INSERT INTO %s VALUES(?,?,?,?,?,?)", tableName);

        try {
            pstmt = conn.prepareStatement(query);

            pstmt.setInt(1, ep.getT());
            pstmt.setString(2, ep.getBodyName());
            pstmt.setDouble(3, ep.getRx());
            pstmt.setDouble(4, ep.getRy());
            pstmt.setDouble(5, ep.getVx());
            pstmt.setDouble(6, ep.getVy());

            pstmt.executeUpdate();
        }
        catch (SQLException e) {
            System.out.println(e.getMessage());
            throw new RuntimeException("Can't save data to database");
        }
    }

    public static void saveEphemeris(String tableName, Planet pl, int daysPassed) {
        //doesn't need to add the data if it already exists
        if(rowExists(tableName, pl.getName(), daysPassed)) {
            return;
        }

        String query = String.format("INSERT INTO %s VALUES(?,?,?,?,?,?)", tableName);

        try {
            pstmt = conn.prepareStatement(query);

            pstmt.setInt(1, daysPassed);
            pstmt.setString(2, pl.getName());
            pstmt.setDouble(3, pl.getLoc().x);
            pstmt.setDouble(4, pl.getLoc().y);
            pstmt.setDouble(5, pl.getVel().x);
            pstmt.setDouble(6, pl.getVel().y);

            pstmt.executeUpdate();
        }
        catch (SQLException e) {
            System.out.println(e.getMessage());
            throw new RuntimeException("Can't save data to database");
        }
    }

    public static boolean rowExists(String tableName, String nameOfPlanet, int daysPassed) {
        String query = String.format("SELECT (count(*) > 0) FROM %s WHERE body_name='%s' AND t=%s",
                tableName,
                nameOfPlanet,
                daysPassed
        );
        try {
            pstmt = conn.prepareStatement(query);

            try {
                // Only expecting a single result
                resSet = pstmt.executeQuery();
                if (resSet.next()) {
                    return resSet.getBoolean(1);
                }
                else {
                    throw new RuntimeException("Can't find the row");
                }
            }
            catch(Exception e) {
                System.out.println(e.getMessage());
                throw new RuntimeException("Can't execute query");
            }
        }
        catch (SQLException e) {
            System.out.println(e.getMessage());
            throw new RuntimeException("Can't check if row exists");
        }
    }

    public static Vector[] getVectorsByNameAndTime(String tableName, String nameOfPlanet, int daysPassed) {
        Vector[] res = new Vector[2];
        //System.out.println(nameOfPlanet + " " + daysPassed);

        String query = String.format("SELECT rx,ry,vx,vy FROM %s WHERE body_name = '%s' AND t = %s",
                tableName,
                nameOfPlanet,
                daysPassed
        );

        try {
            pstmt = conn.prepareStatement(query);

            try {
                resSet = pstmt.executeQuery();
                Vector r = new Vector(
                        resSet.getDouble(1),
                        resSet.getDouble(2)
                );
                Vector v = new Vector(
                        resSet.getDouble(3),
                        resSet.getDouble(4)
                );

                res[0] = r;
                res[1] = v;
            }
            catch(Exception e) {
                //System.out.println("CANT: " + nameOfPlanet + " " + daysPassed);
                System.out.println(e.getMessage());
                throw new RuntimeException("Can't execute query");
            }
        }
        catch (SQLException e) {
            System.out.println(e.getMessage());
            throw new RuntimeException("Can't get vectors from database");
        }

        return res;
    }

    public static void closeDB() throws SQLException {
        conn.close();
        pstmt.close();
        resSet.close();
    }
}
