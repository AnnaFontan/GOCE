import java.lang.Math;

public class OP {
    private OP() {
        // restrict instantiation
    }

    public static double norm(double[] array) {
        if (array != null) {
            double sum = 0.0;
            for (int i = 0; i < array.length; i++) {
                sum += Math.pow(array[i], 2);
            }
            return Math.pow(sum, 0.5);
        }
        else { return 0.0; }
    }

    public static double[] array_coeff_product(double coeff, double[] array) {
        double[] prod = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            prod[i] = array[i] * coeff;
        }
        return prod;
    }

    public static double array_scalar_product(double[] array1, double[] array2) {
        double prod = 0;
        for (int i = 0; i < array1.length; i++) {
            prod += array1[i] * array2[i];
        }
        return prod;
    }

    public static double[] array_cross_product(double[] array1, double[] array2) {
        // Check that the two array have the same length
        double[] prod = new double[array1.length];
        prod[0] = array1[1]*array2[2] - array1[2]*array2[1];
        prod[1] = array1[2]*array2[0] - array1[0]*array2[2];
        prod[2] = array1[0]*array2[1] - array1[1]*array2[0];
        return prod;
    }

    public static double[] array_sum(double[] array1, double[] array2) {
        // Check that the two array have the same length
        double[] sum = new double[array1.length];
        for (int i = 0; i < array1.length; i++) {
            sum[i] = array1[i] + array2[i];
        }
        return sum;
    }

    public static double[] array_sub(double[] array1, double[] array2) {
        // Check that the two array have the same length
        double[] sum = new double[array1.length];
        for (int i = 0; i < array1.length; i++) {
            sum[i] = array1[i] - array2[i];
        }
        return sum;
    }

    public static double[][] matrix_product(double[][] A, double[][] B) {
        double[][] C = new double[A.length][B[0].length];
        double sum;

        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < B[0].length; j++) {
                sum = 0;
                for (int r = 0; r < A[0].length; r++) {
                    sum += A[i][r] * B[r][j];
                }
                C[i][j] = sum;
            }
        }
        return C;
    }

    public static double[][] matrix_transpose(double[][] A) {
        double[][] T = new double[A[0].length][A.length];

        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                T[j][i] = A[i][j];
            }
        }
        return T;
    }

    public static double[] matrix_array_product(double[][] A, double[] array) {
        double[] b = new double[array.length];
        double sum;

        for (int i = 0; i < A.length; i++) {
            sum = 0;
            for (int j = 0; j < A[0].length; j++) {
                sum += A[i][j] * array[j];
            }
            b[i] = sum;
        }
        return b;
    }
}