import java.util.*;
import org.apache.commons.math3.distribution.*;

/**
 * Let me use this for the main test harness.
 * 
 * @author Daniel Seita
 *
 */
public class MainTestHarness {
    
    private static Random rand = new Random();

    /**
     * This will run a bunch of experiments and print to standard output, so it's probably best to
     * save all this info in random files.
     * 
     * @param args Nothing for now
     */
    public static void main(String[] args) {
        System.out.println("Now testing the single Page's test with one state test and Gaussians ...");
        testSinglePageOneState();
        System.out.println("... done.");
    }
    
    /**
     * The purpose of this is to test out a single sensor when the data is generated from a single source,
     * from two different Gaussian distributions, where the threshold that the Gaussian distribution changes
     * is set at the 500^th-indexed iteration throughout all trials. We only end trials after we see a score
     * that exceeds the threshold. False positives are counted when Page's test exceeds the gamma threshold
     * before the 500^th-indexed iteration. Right now this should be done, though we might tweak some values
     * to make the resulting graphs/plots more readable.
     */
    public static void testSinglePageOneState() {
        NormalDistribution n1 = new NormalDistribution(0.0, 1.0);
        NormalDistribution n2 = new NormalDistribution(0.0, 1.032);
        for (double threshold = 1.0; threshold <= 10; threshold++) {
            int totalFalseAlarms = 0;   // Total times Page's test exceeds thresh when it's benign
            int totalSamplesFalse = 0;  // Total samples drawn before malignant distribution begins
            int totalSamples = 0;       // Total samples drawn overall across all trials here
            double averageDetectionDelay = -1.0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= 10000; trial++) {
                // Now let's go through the iterations to perform the test
                boolean pageTestDone = false;
                int threshForChange = 500; // Let's do this for simplicity.
                int iteration = 0;
                double previousScore = 0.0;
                while (!pageTestDone) {
                    double observation = 0.0;
                    if (iteration < threshForChange) {
                        totalSamplesFalse++;
                        observation = n1.sample();
                    } else {
                        observation = n2.sample();
                    }
                    double currentScore = Math.max(0.0, 
                            previousScore + Math.log(n2.density(observation) / n1.density(observation)));
                    previousScore = currentScore;
                    if (currentScore > threshold) {
                        if (iteration < threshForChange) {
                            totalFalseAlarms++;
                        } else if (iteration >= threshForChange) {
                            pageTestDone = true;
                            int thisDetectionDelay = iteration - 500;
                            averageDetectionDelay = ((trial-1)*averageDetectionDelay + thisDetectionDelay) / ((double) trial);
                        }
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            System.out.println("False alarms = " + totalFalseAlarms);
            System.out.println("Total samples (before thresh) = " + totalSamplesFalse);
            System.out.println("Total samples = " + totalSamples);
            double ratio = ((double) totalFalseAlarms) / totalSamplesFalse;
            System.out.println("Ratio = " + ratio);
            System.out.println("Average detection delay = " + averageDetectionDelay);
        }
    }

    /**
     * Returns a pseudo-random number between min and max, inclusive. The difference between 
     * min and max can be at most <code>Integer.MAX_VALUE - 1</code>.
     *
     * @param min Minimum value
     * @param max Maximum value. Must be greater than min.
     * @return Integer between min and max, inclusive.
     * @see java.util.Random#nextInt(int)
     */
    public static int randInt(int min, int max) {
        int randomNum = rand.nextInt((max - min) + 1) + min;
        return randomNum;
    }
}
