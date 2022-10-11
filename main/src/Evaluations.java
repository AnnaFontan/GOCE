import java.lang.Math;
public class Evaluations {
    private Evaluations() {
        // restrict instantiation
    }

    public static double ion_thruster(double m_dot) {
        // Function used to evaluate the thrust and to analyse the ion thruster system
        double T = m_dot * Math.sqrt(2*Constants.IonT.e*Constants.IonT.DV/Constants.IonT.mi) * 1e-3; // [kN] thrust
        return T;
    }


    public static double[] accelerator(double x_a, double v_a, double Vout, double D, double T, int j) {
        // Function to solve the system related to the accelerometer

        double dVout = 2*v_a*Constants.acc.Vbias*Constants.acc.epsilon*Constants.acc.Aa/Constants.acc.Cf*
                (Math.pow(Constants.acc.g, 2) + Math.pow(x_a, 2))/Math.pow((Math.pow(Constants.acc.g, 2) - Math.pow(x_a, 2)), 2); // [V/s] dot(Vout)

        double Vc = Constants.acc.Kpa*Vout + Constants.acc.Kda*dVout; // [V]: PD controller
        double Vx = x_a/Constants.acc.g*Constants.acc.Vbias; // [V]

        double DV_1 = Constants.acc.Vbias - Vc - 0.5*Vx; // [V] first difference voltages
        double DV_2 = Constants.acc.Vbias + Vc + 0.5*Vx; // [V] second difference voltages

        double F_ext = (T - D)*1e3*Constants.acc.m/Constants.GOCE.M; // [N] external force
        double F_el1 = 0.5 * Constants.acc.epsilon*Constants.acc.Aa*Math.pow(DV_1/(Constants.acc.g - x_a), 2); // [N] first electrostatic force
        double F_el2 = 0.5 * Constants.acc.epsilon*Constants.acc.Aa*Math.pow(DV_2/(Constants.acc.g + x_a), 2); // [N] second electrostatic force

        double[] dY = new double[8];
        if (j == 0) { // NO LINEARISATION
            dY[0] = v_a; // dY(1) = dot(x_a)
            dY[1] = -(F_ext - F_el1 + F_el2) / Constants.acc.m; // dY(2) = dot(v_a)
            dY[2] = dVout; // dY(3) = dot(Vout)
            dY[3] = (Constants.valve.Kpv * dVout + Constants.valve.Kiv * Vout); // dY(4) = dot(I):PI controller
            // out_acc = 0;
        }
        else { // linearisation
            // dY = 0; // don't need these values in the linearisation case
            dY[4] = F_ext; // out_acc[0]
            dY[5] = F_el1; // out_acc[1]
            dY[6] = F_el2; // out_acc[2]
            dY[7] = dVout; // out_acc[3]
        }

        return dY;
    }



    public static double[] flow_control_valve(double IonT_P2, double I, double x_v, double v_v, int j) {
        // Function to solve the flow control valve system

        if (j == 0) { // NO LINEARISATION
            /* Valve boundary conditions:
               for x_v = 0: the valve is fully open --> A_d = valve.A_0
               for x_v = 2*valve.r0: the valve is fully closed --> A_d = 0
             */
            if (x_v < 0) { x_v = 0; }
            if (x_v <= 0 && v_v < 0) { x_v = 0; v_v = 0; }
            if (x_v > 2*Constants.valve.r0) { x_v = 2*Constants.valve.r0; }
            if (x_v >= 2*Constants.valve.r0 && v_v > 0) { x_v = 2*Constants.valve.r0; v_v = 0; }
        }

        double z = (2*Constants.valve.r0 - x_v)/(2*Constants.valve.r0); // [-] valve command, goes from 0 to 1
        double alpha = 2*Math.acos(1 - 2*z); // [rad] central angle measuring valve opening
        double A_d = Math.pow(Constants.valve.r0, 2)/2 * (alpha - Math.sin(alpha)); // [m^2] orifice area
        double m_dot = IonT_P2*A_d/Math.pow((Constants.IonT.k + 1)/2, (Constants.IonT.k + 1)/(2*(Constants.IonT.k - 1))) *
                Math.sqrt(Constants.IonT.k/(Constants.IonT.Rm*Constants.IonT.T2)); // [kg/s] mass flow rate

        double[] dY = new double[4];
        if (j == 0) {
            // State derivatives
            dY[0] = v_v; // dY(11) = dot(x_v)
            dY[1] = 1/Constants.valve.m_fcv * (Constants.valve.Kfcv*(2*Constants.valve.r0 - x_v) - Constants.valve.c*v_v - Constants.valve.Ki*I); // dY(12) = dot(v_v)
        }

        dY[2] = A_d;
        dY[3] = m_dot;

        return dY;
    }


    public static double[] orbital_mechanics (double[] kepl, double T, double[] D_v_in, int j) {
        // Function to solve the Orbital Mechanic system

        double h = kepl[0]; // keplerian parameters
        double e = kepl[1];
        double in = kepl[2];
        double OM = kepl[3];
        double om = kepl[4];
        double th = kepl[5];

        double u = th + om; // [rad] argument of latitude
        double sma = Math.pow(h, 2)/(Constants.Earth.mu*(1 - Math.pow(e, 2))); // [km] semi-major axis
        [RR, VV] = kep2car(sma, e, in, OM, om, th, Constants.Earth.mu); // radius and velocity vectors

        /* Earth's radius (2D orbit) has been evaluated in cartesian coordinate using kep2car.
           Then the ellipse has been translated into the center of the Earth.
         */
        double E = th + om;
        double thE = 2*Math.atan(Math.sqrt((1 + Constants.Earth.e)/(1 - Constants.Earth.e))*Math.tan(E/2)); // [rad] angle of GOCE projection on Earth's surface
        [RR_Earth, ~] = kep2car(Constants.Earth.sma, Constants.Earth.e, Constants.orbit.i0, Constants.orbit.OM0, Constants.orbit.om0, thE, Constants.Earth.mu); // Earth's radius
        RR_Earth[0] = RR_Earth[0] + Constants.Earth.sma*Constants.Earth.e; // translation into Earth's center

        double Rn = Math.pow(h, 2)/(Constants.Earth.mu*(1 + e*Math.cos(th))); // [km] radius norm
        double Rx = RR[0]; // [km] radius in x-direction
        double Ry = RR[1]; // [km] radius in y-direction
        double Rz = RR[2]; // [km] radius in z-direction

        // Perturbing gravitation acceleration: J2
        double[] array = { Rx/Rn*(5*Math.pow(Rz/Rn, 2) - 1), Ry/Rn*(5*Math.pow(Rz/Rn, 2) - 1), Rz/Rn*(5*Math.pow(Rz/Rn, 2) - 3) };
        double p_J2_coeff = 1.5*(Constants.Earth.J2*Constants.Earth.mu*Math.pow(Constants.Earth.sma, 2))/(Math.pow(Rn, 4));
        double[] p_J2 = OP.array_coeff_product(p_J2_coeff, array);

        // The values of default are the ones used for the linearisation
        double[] p_drag = OP.array_coeff_product(1/Constants.GOCE.M, D_v_in); // [kN]drag force (vector)
        double H = Constants.orbit.alt0; // fixed to its initial value
        double rho = 0;
        double D = OP.norm(D_v_in); // in the linearisation the Drag is fixed to its initial value
        double[] D_v = {0, 0, 0};

        if (j == 0) { // NO LINEARISATION
            double[] v_rel = OP.array_sub(VV, OP.array_cross_product(Constants.Earth.wE, RR)); // [km / s]relative velocity
            H = OP.norm(OP.array_sub(RR, RR_Earth)); // [km]altitude
            rho = Coversions.density(H); // [kg / km ^ 3]

            // Perturbing gravitation acceleration: drag
            p_drag = OP.array_coeff_product(-0.5*rho*OP.norm(v_rel)/ Constants.GOCE.beta, v_rel); // [km / s^2] acceleration due to drag
            D_v = OP.array_coeff_product(Constants.GOCE.M, p_drag); // [kN]drag force (vector)
            D = OP.norm(D_v); // [kN] drag(scalar)
        }

        double[] p_thrust = OP.array_coeff_product(-T/Constants.GOCE.M/OP.norm(p_drag), p_drag);

        // Unitary vectors
        double[] r = { - Math.sin(OM)*Math.cos(in)*Math.sin(u)+Math.cos(OM)*Math.cos(u), Math.cos(OM)*Math.cos(in)*Math.sin(u)+Math.sin(OM)*Math.cos(u), Math.sin(in)*Math.sin(u) }; // (3x1)
        double[] s = { - Math.sin(OM)*Math.cos(in)*Math.cos(u)-Math.cos(OM)*Math.sin(u), Math.cos(OM)*Math.cos(in)*Math.cos(u)-Math.sin(OM)*Math.sin(u), Math.sin(in)*Math.cos(u) }; // (3x1)
        double[] w = { Math.sin(OM)*Math.sin(in), - Math.cos(OM)*Math.sin(in), Math.cos(in) }; // (3x1)

        // Components of the perturbing accelerations in the LVLH frame
        double[] vect = OP.array_sum(OP.array_sum(p_J2, p_drag), p_thrust);
        double pr = OP.array_scalar_product(vect, r); // (1x1)
        double ps = OP.array_scalar_product(vect, s); // (1x1)
        double pw = OP.array_scalar_product(vect, w); // (1x1)

        // The values of default are the ones of the linearisation case
        double[] dY = {0, 0, 0, 0, 0}; // don't need these values in the linearisation case
        double[] out_lin = {Rn, pr, ps, pw}; // output parameters
        if (j == 0) {
        // Orbital Mechanics: first derivatives
            dY[0] = Rn*ps; // dY(1) = dot(h)
            dY[1] = h/Constants.Earth.mu*Math.sin(th)*pr + 1/(Constants.Earth.mu*h)*((Math.pow(h, 2) + Constants.Earth.mu*Rn)*Math.cos(th) + Constants.Earth.mu*Rn*e)*ps; // dY(2) = dot(e)
            dY[2] = Rn/h*Math.cos(u)*pw; // dY(3) = dot(i)
            dY[3] = Rn*Math.sin(u)/(h*Math.sin(in))*pw; // dY(4) = dot(OM)
            dY[4] = - 1/(e*h) * (Math.pow(h, 2)/Constants.Earth.mu*Math.cos(th)*pr - (Rn + Math.pow(h, 2)/Constants.Earth.mu)*Math.sin(th)*ps) - Rn*Math.sin(u)/(h*Math.tan(in))*pw; // dY(5) = dot(om)
            dY[5] = h/Math.pow(Rn, 2) + 1/(e*h)*(Math.pow(h, 2)/Constants.Earth.mu*Math.cos(th)*pr - (Math.pow(h, 2)/Constants.Earth.mu + Rn)*Math.sin(th)*ps); // dY(6) = dot(th)
            out_lin = OP.array_coeff_product(0, out_lin); // this vector is used only in the linearisation case
        }

        // Output of the function
        dY[6] = D;
        dY[7] = rho;
        dY[8] = H;
        dY[9] = sma;

        dY[10] = D_v[0];
        dY[11] = D_v[1];
        dY[12] = D_v[2];

        dY[13] = out_lin[0];
        dY[14] = out_lin[1];
        dY[15] = out_lin[2];
        dY[16] = out_lin[3];

        dY[17] = RR[0];
        dY[18] = RR[1];
        dY[19] = RR[2];

        dY[20] = RR_Earth[0];
        dY[21] = RR_Earth[1];
        dY[22] = RR_Earth[2];
    }
}
