// Function to define the data of the project
import java.lang.Math;
public final class Constants {
    public static final double Ru = 8.314462; // [J/(K mol)] universal gas constant

    // GOCE
    public final class GOCE {
        public static final double M = 300.0; // [kg] GOCE mass
        public static final double A = 1.1; // [m^2] GOCE frontal section
        public static final double beta = 300*1e6; // [kg/km^2] GOCE ballistic coefficient
    }
    // Orbital Mechanic data
    public final class Earth {
        public static final double sma = 6378.16; // [km] semi-major axis: shape of Earth
        public static final double b = 6356.778; // [km] semi-minor axis: shape of Earth
        public static final double mu = 398600; // [km^3/s^2] Earth gravitational constant
        public static final double[] wE = {0, 0, (15.04 * Math.PI / 180) / 3600}; // [rad/s] Earth rotation ratio
        public static final double J2 = 0.00108263; // Earth oblateness effect
        public static final double e = Math.sqrt(1 - Math.pow(b / sma, 2)); // [-] eccentricity of the Earth
    }

    // Initial conditions
    public final class orbit {
        public static final double alt0 = 254.9; // [km] altitude of the orbit
        public static final double e0 = 0.0045; // [-] eccentricity of the orbit
        public static final double i0 = Math.PI/2; // [rad] inclination of the orbit
        public static final double th0 = 0; // [rad] true anomaly
        public static final double om0 = 0; // [rad] anomaly of the perigee
        public static final double OM0 = 0; // [rad] right ascension of the ascending node
        public static final double rp = alt0 + Earth.sma; // [km] radius of the pericenter
        public static final double sma0 = rp/(1 - e0); // [km] semi-major axis
        public static final double T = 2*Math.PI*Math.sqrt(Math.pow(sma0, 3)/Earth.mu); // [s] period of the orbit
        public static final double h0 = Math.sqrt(Earth.mu*sma0*(1 - Math.pow(e0, 2))); // [m^2/s] specific angular momentum
    }
    // Accelerometer data
    public final class acc {
        public static final double epsilon = 8.85*1e-12; // [F/m] accelerometer permittivity
        public static final double Aa = 1.6*1e-3; // [m^2] accelerometer seismic mass section
        public static final double m = 0.32; // [kg] accelerometer seismic mass
        public static final double g = 5*1e-4; // [m] electrodes - mass gap
        public static final double Vbias = 10; // [V] accelerometer bias voltage
        public static final double Cf = 2*1e-12; // [F] capacitance
        public static final double Kpa = 1e6; // [-] PD controller proportional gain --> OPTIMIZABLE
        public static final double Kda = 5*1e4; // [1/s] PD controller derivative gain --> OPTIMIZABLE
        public static final double L = 1e-3; // [H] inductance of the solenoid
    }
    // Control valve data
    public final class valve {
        public static final double A0 = 1e-5; // [m^2] valve orifice area
        public static final double r0 = Math.sqrt(A0/Math.PI); // [m] valve orifice radius
        public static final double m_fcv = 2*1e-1; // [kg] valve spool mass --> OPTIMIZABLE
        public static final double Ki = 0.2; // [-] proportionality coefficient current-spool --> OPTIMIZABLE
        public static final double Kpv = 0.1; // [-] PI controller proportional gain --> OPTIMIZABLE
        public static final double Kiv = 3; // [s] PI controller integral gain --> OPTIMIZABLE
        public static final double Kfcv = 7*1e3; // [N/m] valve spring coefficient
        public static final double c = 30; // [N s/m] valve friction coefficient
    }
    // Ion thruster data
    public final class IonT {
        public static final double T2 = 240; // [K] Xenon working temperature
        public static final double P2_in = 2*1e5; // [Pa] Xenon working pressure
        public static final double k = 1.66; // [-] Xenon specific heat ratio
        public static final double mi = 2.188*1e-25; // [kg] Xenon ion mass
        public static final double e = 1.6*1e-19; // [C] electron Charge
        public static final double DV = 2*1e3; // [V] acceleration grid voltage
        public static final double Mm = 131.293*1e-3; // Xenon molar mass [kg/mol]
        public static final double Rm = Ru/Mm; // [J/(kg K)] specific gas constant
    }
    // Failures
    public final class failure {
        public static final double I = 2*3600; // [s] failure of the current happens after two hours
        public static final double T = 3600; // [s] failure of the thrust happens after one hour
        public static final double PXe_low = 12*3600; // [s] the PXe reaches its lowest value after 12hrs
    }
}