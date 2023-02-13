import java.util.*;
import java.io.*;
import java.text.*;

// X_0, Y_0, X_1, Y_1
public class PointScatteringRep1 extends FitnessFunction {
    private static final int BITS = 21;

    public PointScatteringRep1() throws IOException {
        name = "Point Scattering Problem Representation #1";
    }

    private double dist(double[] p1, double[] p2) {
        return Math.sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]));
    }

    private void mapToUnitCircle(double[] p) {
        for (int i = 0; i < p.length; i++) {
            p[i] -= (1 << (BITS - 1));
            p[i] /= 1e6;

            p[i] = Math.min(p[i], 1.0);
            p[i] = Math.max(p[i], -1.0);
        }
    }

    public void doRawFitness(Chromo chromo) {
        double maxDist = -1;
        for (int i = 0; i < Parameters.numGenes - 2; i += 2) {
            double[] p_i = new double[]{chromo.getPosIntGeneValue(i), chromo.getPosIntGeneValue(i + 1)};
            mapToUnitCircle(p_i);

            double minDist = 1e9;
            for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                double[] p_j = new double[]{chromo.getPosIntGeneValue(j), chromo.getPosIntGeneValue(j + 1)};
                mapToUnitCircle(p_j);

                // System.out.println("DEBUG: p_i: " + p_i[0] + " " + p_i[1] + " p_j: " + p_j[0] + " " + p_j[1] + " dist: " + dist(p_i, p_j));

                minDist = Math.min(minDist, dist(p_i, p_j));
            }

            maxDist = Math.max(maxDist, minDist);
        }

        chromo.rawFitness = maxDist;
    }
}