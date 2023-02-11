import java.util.*;
import java.io.*;
import java.text.*;
import java.awt.Point;

// X_0, Y_0, X_1, Y_1
public class PointScatteringRep1 extends FitnessFunction {
    public PointScattering() throws IOException {
        name = "Point Scattering Problem";
    }

    private double dist(Point p1, Point p2) {
        return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
    }

    public void doRawFitness(Chromo chromo) {
        double maxDist = -1;
        for (int i = 0; i < Parameters.numGenes; i += 2) {
            Point p_i = new Point(chromo.getPosIntGeneValue(i) / 1e6, chromo.getPosIntGeneValue(i + 1) / 1e6);

            double minDist = 1e9;
            for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                Point p_j = new Point(chromo.getPosIntGeneValue(j) / 1e6, chromo.getPosIntGeneValue(j + 1) / 1e6);

                minDist = Math.min(minDist, dist(p_i, p_j));
            }

            maxDist = max(maxDist, minDist);
        }

        chromo.rawFitness = maxDist;
    }
}