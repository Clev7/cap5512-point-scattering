import java.util.*;
import java.io.*;
import java.text.*;

// X_0, X_1, Y_0, Y_1
public class PointScatteringRep2 extends FitnessFunction {
    private static final int BITS = 21;

    public PointScatteringRep2() throws IOException {
        name = "Point Scattering Problem Representation #2";
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
        int n = Parameters.numGenes / 2;
        for (int i = 0; i < n - 1; i++) {
            double[] p_i = new double[]{chromo.getPosIntGeneValue(i), chromo.getPosIntGeneValue(i + n)};
            mapToUnitCircle(p_i);

            double minDist = 1e9;
            for (int j = i + 1; j < n; j++) {
                double[] p_j = new double[]{chromo.getPosIntGeneValue(j), chromo.getPosIntGeneValue(j + n)};
                mapToUnitCircle(p_j);

                minDist = Math.min(minDist, dist(p_i, p_j));
            }

            maxDist = Math.max(maxDist, minDist);
        }

        chromo.rawFitness = maxDist;
    }
}
