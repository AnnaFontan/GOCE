import java.io.*;
import org.apache.commons.math3.ode.*;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.sampling.StepInterpolator;
public class Main {
    public static void main(String[] args) {

        double[] y0 = {Constants.orbit.h0, Constants.orbit.e0, Constants.orbit.i0, Constants.orbit.OM0, Constants.orbit.om0,
                Constants.orbit.th0, 0, 0, 0, 0, 2*Constants.valve.r0, 0}; // vector of initial conditions
        double[] tspan0 = {0, 6*Constants.orbit.T}; // time vector [s]

        // 1st integration
        String failure = ""; // in this case, no failures considered
        double start = System.nanoTime(); // tic (Matlab)

        FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1e-14, 1e-5, 1e-14, 1e-13);
        FirstOrderDifferentialEquations fun = new ode_GOCE();

        dp853.integrate(fun, tspan0[0], y0, tspan0[1], y0);

        System.out.println(y0[0]);

        double end = System.nanoTime();
        double t_firstI = (end - start)*1e-9; // [s] evaluation of the time needed to compute the integration

        System.out.println("1. First integration completed: t = [" + tspan0[0]/3600 + ", " + tspan0[1]/3600 + "] hrs");
        System.out.println("\t - duration = " + t_firstI + " s \n");

    }
}