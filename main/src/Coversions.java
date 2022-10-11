import java.lang.Math;

public class Coversions {
    private Coversions() {
        // restrict instantiation
    }

    public static double density(double H) {
        /* Function that evaluated the density, given the altitude. See reference [1] in the PDF.
           Model: DTM-2013 thermosphere model. NB: the function works only for altitudes between 100 to 900 km
         */

        // Coefficient from the data fitting
        double A = -4.467926544726752;
        double B = -0.08731069078582779;
        double C = 0.00038339937759263866;
        double D = -0.00000102392405619931238;
        double E = 1.4872119104109197e-9;
        double F = -1.0821365008660183e-12;
        double G = 3.115172734059814e-16;

        double rho = Math.pow(10, A + B*H + C*Math.pow(H, 2) + D*Math.pow(H, 3) + E*Math.pow(H, 4) + F*Math.pow(H, 5) + G*Math.pow(H, 6)) * 1e12; // interpolation of the density [kg/km^3]

        return rho;
    }


    public static double[] kep2car(double a, double e, double i, double OM, double om, double th, double mu) {
        // Function used in order to obtain the keplerian coordinates from the cartesian ones

        int j;

        double p = a * (1 - Math.pow(e,2)); // [km]semi - latus rectus
        double r = p / (1 + e * Math.cos(th)); // [km]radius

        // Perifocal frame
        double[] rPF_direction = {Math.cos(th), Math.sin(th), 0};
        double[] vPF_direction = {-Math.sin(th), e + Math.cos(th), 0};
        double[] rPF = OP.array_coeff_product(r, rPF_direction); // [km] radius vector
        double[] vPF = OP.array_coeff_product(Math.sqrt(mu / p), vPF_direction); // [km / s] velocity vector

        // Rotation matrices
        double[][] R3OM = {{Math.cos(OM), Math.sin(OM), 0}, {-Math.sin(OM), Math.cos(OM), 0}, {0, 0, 1}};
        double[][] R1i = {{1, 0, 0}, {0, Math.cos(i), Math.sin(i)}, {0, -Math.sin(i), Math.cos(i)}};
        double[][] R3om = {{Math.cos(om), Math.sin(om), 0}, {-Math.sin(om), Math.cos(om), 0}, {0, 0, 1}};

            GE2PF = R3om * R1i * R3OM; // cartesian ---> perifocal
            PF2GE = GE2PF.'; // perifocal ---> cartesian

            // Cartesian frame
            rGE = PF2GE * rPF; // [km]radius vector
            vGE = PF2GE * vPF; // [km / s]velocity vector

            rr = rGE; // [km]radius vector
            vv = vGE; // [km / s]velocity vector
        }
    }
}