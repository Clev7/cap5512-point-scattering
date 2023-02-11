import java.util.*;
import java.io.*;
import java.text.*;

// X_0, Y_0, X_1, Y_1
public class PointScatteringRep1 extends FitnessFunction {
    public PointScatteringRep1() throws IOException {
        name = "Point Scattering Problem Representation #1";
    }

    private double dist(double[] p1, double[] p2) {
        return Math.sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]));
    }

    public void doRawFitness(Chromo chromo) {
        double maxDist = -1;
        for (int i = 0; i < Parameters.numGenes; i += 2) {
            double[] p_i = new double[]{chromo.getPosIntGeneValue(i) / 1e6, chromo.getPosIntGeneValue(i + 1) / 1e6};

            double minDist = 1e9;
            for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                double[] p_j = new double[]{chromo.getPosIntGeneValue(j) / 1e6, chromo.getPosIntGeneValue(j + 1) / 1e6};

                minDist = Math.min(minDist, dist(p_i, p_j));
            }

            maxDist = Math.max(maxDist, minDist);
        }

        chromo.rawFitness = maxDist;
    }
}