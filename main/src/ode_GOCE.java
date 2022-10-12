import java.lang.Math;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
public class ode_GOCE implements FirstOrderDifferentialEquations {

    public int getDimension() {
        return 12;
    }

    public void computeDerivatives(double t, double[] y, double[] dY) {
        // Function of first order ODEs to solve the whole system

        String failure = "";
        int linearisation = 0; // no linearisation considered

        // Orbital Mechanics
        double h = y[0]; // [km^2/s] specific angular momentum
        double e = y[1]; // [-] eccentricity
        double in = y[2]; // [rad] inclination
        double OM = y[3]; // [rad] right ascension of the ascending node
        double om = y[4]; // [rad] argument of the perigee
        double th = y[5]; // [rad] true anomaly
        // Accelerometer + Control Valve
        double x_a = y[6]; // [m] accelerometer position
        double v_a = y[7]; // [m/s] accelerometer velocity
        double Vout = y[8]; // [V] output voltage
        double I = y[9]; // [A] current
        double x_v = y[10]; // [m] valve position
        double v_v = y[11]; // [m/s] valve velocity

        // Valve equations ---------------------------------------------------------------------------------------------
        double IonT_P2;
        String s1 = "PXe";
        if (s1.equals(failure)) { // if failure associated with the pressure P_Xe
            double m = -0.9 * Constants.IonT.P2_in / Constants.failure.PXe_low;
            IonT_P2 = Constants.IonT.P2_in + m * t; // the pressure decreases linearly till the 10% of its starting value
        } else {
            IonT_P2 = Constants.IonT.P2_in;
        }
        double[] dY1 = Evaluations.flow_control_valve(IonT_P2, I, x_v, v_v, linearisation);
        double[] dY_Valve = new double[2];
        System.arraycopy(dY1, 0, dY_Valve, 0, dY_Valve.length);
        double A_d = dY1[2];
        double m_dot = dY1[3];

        // Ion Thruster equations --------------------------------------------------------------------------------------
        double T = Evaluations.ion_thruster(m_dot);
        String s2 = "T";
        if (t >= Constants.failure.T && s2.equals(failure)) {// if failure associated with the T
            T = T * Math.exp(-10 * (t - Constants.failure.T) / Constants.failure.T); // temperature decreasing exponentially
        }

        // Orbital Mechanics equations ---------------------------------------------------------------------------------
        double[] kepl = {h, e, in, OM, om, th};
        double[] D_V_in = {0, 0, 0}; // for linearisation = 0
        double[] dY2 = Evaluations.orbital_mechanics(kepl, T, D_V_in, linearisation);
        double[] dY_OrbMech = new double[6];
        System.arraycopy(dY2, 0, dY_OrbMech, 0, dY_OrbMech.length);
        double D = dY2[6];
        double rho = dY2[7];
        double H = dY2[8];
        double sma = dY2[9];
        double[] D_v = {dY2[10], dY2[11], dY2[12]};
        double[] RR = {dY2[17], dY2[18], dY2[19]};
        double[] RR_Earth = {dY2[20], dY2[21], dY2[22]};

        // Accelerometer equations -------------------------------------------------------------------------------------
        double[] dY3 = Evaluations.accelerator(x_a, v_a, Vout, D, T, linearisation);
        double[] dY_Acc = new double[4];
        System.arraycopy(dY3, 0, dY_Acc, 0, dY_Acc.length);

        // State vector: differential equations ------------------------------------------------------------------------
        String s3 = "I";
        System.arraycopy(dY_OrbMech, 0, dY, 0, dY_OrbMech.length); // differential equations related to the Orbital Mechanic system [6x1]
        System.arraycopy(dY_Acc, 0, dY, dY_OrbMech.length, dY_Acc.length); // differential equations related to the system of the accelerometer [4x1]
        System.arraycopy(dY_Valve, 0, dY, dY_OrbMech.length + dY_Acc.length, dY_Valve.length); // differential equations related to the system of the valve [2x1]

        if (t >= Constants.failure.I && s3.equals(failure)) {// if failure associated with the current I
            dY[9] = -I * 1e-3; // exponential decrease of the current
        }

        // Output parameters
        double[] output = new double[28];
        System.arraycopy(dY, 0, output, 0, dY.length);
        output[12] = sma;
        output[13] = H;
        output[14] = rho;
        output[15] = D;
        output[16] = T;
        output[17] = A_d;
        output[18] = m_dot;
        System.arraycopy(D_v, 0, output, 19, 3);
        System.arraycopy(RR, 0, output, 22, 3);
        System.arraycopy(RR_Earth, 0, output, 25, 3);
    }
}