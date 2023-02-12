import java.util.*;
import java.io.*;
import java.text.*;

// X_0, X_1, Y_0, Y_1
public class PointScatteringRep2 extends FitnessFunction {
    public PointScatteringRep2() throws IOException {
        name = "Point Scattering Problem Representation #2";
    }

    private double dist(double[] p1, double[] p2) {
        return Math.sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]));
    }

    public void doRawFitness(Chromo chromo) {
        double maxDist = -1;
        int n = Parameters.numGenes / 2;
        for (int i = 0; i < n - 1; i++) {
            double[] p_i = new double[]{chromo.getPosIntGeneValue(i) / 1e6, chromo.getPosIntGeneValue(i + n) / 1e6};
            p_i[0] = Math.min(p_i[0], 1.0);
            p_i[1] = Math.min(p_i[1], 1.0);

            double minDist = 1e9;
            for (int j = i + 1; j < n; j++) {
                double[] p_j = new double[]{chromo.getPosIntGeneValue(j) / 1e6, chromo.getPosIntGeneValue(j + n) / 1e6};

                minDist = Math.min(minDist, dist(p_i, p_j));
            }

            maxDist = Math.max(maxDist, minDist);
        }

        chromo.rawFitness = maxDist;
    }
}
