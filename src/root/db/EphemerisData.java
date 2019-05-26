package root.db;

public class EphemerisData {
    private int t;
    private String bodyName;
    private float rx;
    private float ry;
    private float vx;
    private float vy;

    public EphemerisData(int t, String name, float rx, float ry, float vx, float vy) {
        this.t = t;
        this.bodyName = name;
        this.rx = rx;
        this.ry = ry;
        this.vx = vx;
        this.vy = vy;
    }

    public int getT() {
        return t;
    }

    public void setT(int t) {
        this.t = t;
    }

    public String getBodyName() {
        return bodyName;
    }

    public void setBodyName(String bodyName) {
        this.bodyName = bodyName;
    }

    public float getRx() {
        return rx;
    }

    public void setRx(float rx) {
        this.rx = rx;
    }

    public float getRy() {
        return ry;
    }

    public void setRy(float ry) {
        this.ry = ry;
    }

    public float getVx() {
        return vx;
    }

    public void setVx(float vx) {
        this.vx = vx;
    }

    public float getVy() {
        return vy;
    }

    public void setVy(float vy) {
        this.vy = vy;
    }
}
