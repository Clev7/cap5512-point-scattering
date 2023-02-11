import java.util.*;
import java.io.*;
import java.text.*;


// r_0, theta_0, r_1, theta_1, etc...
public class PointScatteringRep3 extends FitnessFunction {
    public PointScatteringRep3() throws IOException {
        name = "Point Scattering Problem Representation #2";
    }

    private double dist(double r1, int theta1, double r2, int theta2) {
        // Law of Cosines
        return Math.sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * Math.cos(theta1 - theta2));
    }

    public void doRawFitness(Chromo chromo) {
        double maxDist = -1;

        for (int i = 0; i < Parameters.numGenes - 2; i += 2) {
            double r_i = chromo.getPosIntGeneValue(i) / 1e6;
            int theta_i = chromo.getPosIntGeneValue(i + 1) % 360;

            double minDist = 1e9;
            for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                double r_j = chromo.getPosIntGeneValue(j) / 1e6;
                int theta_j = chromo.getPosIntGeneValue(j + 1) % 360;

                minDist = Math.min(minDist, dist(r_i, theta_i, r_j, theta_j));
            }

            maxDist = Math.max(maxDist, minDist);
        }

        chromo.rawFitness = maxDist;
    }
}