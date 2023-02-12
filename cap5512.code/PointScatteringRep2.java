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
        for (int i = 0; i < Parameters.numGenes - 2; i += 2) {
            double[] p_i = new double[]{chromo.getPosIntGeneValue(i) / 1e6, chromo.getPosIntGeneValue(i + (Parameters.numGenes/2)) / 1e6};
            p_i[0] = Math.min(p_i[0], 1.0);
            p_i[1] = Math.min(p_i[1], 1.0);

            double minDist = 1e9;
            for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                double[] p_j = new double[]{chromo.getPosIntGeneValue(j) / 1e6, chromo.getPosIntGeneValue(j + 1) / 1e6};

                // System.out.println("DEBUG: p_i: " + p_i[0] + " " + p_i[1] + " p_j: " + p_j[0] + " " + p_j[1] + " dist: " + dist(p_i, p_j));

                minDist = Math.min(minDist, dist(p_i, p_j));
            }

            maxDist = Math.max(maxDist, minDist);
        }

        chromo.rawFitness = maxDist;
    }
}
