package root;

public class Vector {
    public double x, y, z;

    public Vector(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public Vector(double x, double y) {
        this.x = x;
        this.y = y;
        this.z = 0;
    }

    @Override public String toString() {
        return String.format("[%f, %f, %f]", x, y, z);
    }

    public double mag() {
        return Math.sqrt(x*x + y*y + z*z);
    }

    public double magSq() {
        return (x*x + y*y + z*z);
    }

    public void normalize() {
        double m = mag();
        if (m != 0 && m != 1) {
            div(m);
        }
    }

    public Vector cross(Vector v) {
        double crossX = y * v.z - v.y * z;
        double crossY = z * v.x - v.z * x;
        double crossZ = x * v.y - v.x * y;

        return new Vector(crossX, crossY, crossZ);
    }

    public double dot(Vector v) {
        return x*v.x + y*v.y + z*v.z;
    }

    public void add(Vector v) {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    public void add(double n) {
        x += n;
        y += n;
        z += n;
    }

    public void sub(Vector v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    public void mult(double n) {
        x *= n;
        y *= n;
        z *= n;
    }

    public void div(double n) {
        x /= n;
        y /= n;
        z /= n;
    }

    public static Vector add(Vector v1, Vector v2) {
        return new Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }

    public static Vector sub(Vector v1, Vector v2) {
        return new Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }

    public static Vector div(Vector v, double n) {
        return new Vector(v.x/n, v.y/n, v.z/n);
    }

    public static Vector mult(Vector v, double n) {
        return new Vector(v.x*n, v.y*n, v.z*n);
    }



}
