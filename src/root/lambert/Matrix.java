package root.lambert;

public class Matrix {
    int M;
    int N;
    double[][] data;

    public Matrix(int M, int N) {
        this.M = M;
        this.N = N;
        data = new double[M][N];
    }

    public Matrix(double[][] data) {
        M = data.length;
        N = data[0].length;
        this.data = new double[M][N];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                this.data[i][j] = data[i][j];
    }

    public Matrix(double[] data) {
        M = 1;
        N = data.length;
        this.data = new double[1][N];
        for (int i = 0; i < N; i++)
            this.data[0][i] = data[i];
    }

    double toDouble(){
        assert M == 1 & N == 1 :
                "M and N must be equal to 1";
        return data[0][0];
    }

    Matrix transpose() {
        Matrix A = new Matrix(N, M);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[j][i] = this.data[i][j];
        return A;
    }

    Matrix plus(Matrix B) {
        Matrix A = this;
        assert B.M == A.M || B.N == A.N:
                "Can't add matrices: illegal dimensions";
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] + B.data[i][j];
        return C;
    }

    Matrix plus(double num) {
        Matrix A = this;
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] + num;
        return C;
    }

    Matrix minus(Matrix B) {
        Matrix A = this;
        assert B.M == A.M || B.N == A.N:
                "Can't subtract matrices: illegal dimensions";
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] - B.data[i][j];
        return C;
    }

    Matrix divide(double num) {
        Matrix A = this;
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] / num;
        return C;
    }

    Matrix multiply(Matrix B) {
        Matrix A = this;
        assert A.N == B.M:
                "Can't multiply matrices: illegal dimensions";
        Matrix C = new Matrix(A.M, B.N);
        for (int i = 0; i < C.M; i++)
            for (int j = 0; j < C.N; j++)
                for (int k = 0; k < A.N; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }
}
