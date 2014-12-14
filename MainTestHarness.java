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
     * save all this info in random files. When plotting, be sure to use pyplot's y-scale setup if
     * using the false alarm rate, because that depends exponentially on the threshold.
     * 
     * @param args Nothing for now
     */
    public static void main(String[] args) {
        System.out.println("Now testing with a one-state, Gaussian source.");
        int numSensors = 10;
        NormalDistribution n1 = new NormalDistribution(0.0, 1.0);
        NormalDistribution n2 = new NormalDistribution(0.0, 1.032);
        System.out.println("\nNow doing a single Page's test ...");
        testSinglePageOneStateGaussian(n1, n2);
        System.out.println("\nNow doing Page's test in parallel ...");
        testParallelPageOneStateGaussian(n1, n2, numSensors);
        System.out.println("\nNow doing the centralized entity Test ...");
        testIdealCenterOneStateGaussian(n1, n2, numSensors);
        System.out.println("\nNow doing running consensus ...");
        testRunningConsensusOneStateGaussian(n1, n2, numSensors);
        System.out.println("Done.");
    }

    /**
     * This will be perhaps the most interesting part of this code: the Running Consensus scheme.
     * 
     * @param n1
     * @param n2
     * @param M The number of sensors
     */
    public static void testRunningConsensusOneStateGaussian(NormalDistribution n1, NormalDistribution n2, int M) {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();

        System.out.println("To conclude: false alarm rates and detection delays:");
        System.out.println(falseAlarmRates);
        System.out.println(detectionDelays);
    }    

    /**
     * This is implemented similar to the parallel and single sensor tests, except we're allowed to combine
     * knowledge from all nodes together.
     * 
     * @param n1
     * @param n2
     * @param M The number of sensors
     */
    public static void testIdealCenterOneStateGaussian(NormalDistribution n1, NormalDistribution n2, int M) {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.25) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            double averageDetectionDelay = -1.0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= 10000; trial++) {
                boolean pageTestDone = false;
                int threshForChange = 500;
                int iteration = 0;
                double previousScore = 0.0;
                while (!pageTestDone) {
                    double[] observations = new double[M];
                    if (iteration < threshForChange) {
                        totalSamplesFalse++;
                        for (int i = 0; i < M; i++) {
                            observations[i] = n1.sample();
                        }
                    } else {
                        for (int i = 0; i < M; i++) {
                            observations[i] = n2.sample();
                        }
                    }
                    double sumOfLogLikelihoods = 0.0;
                    for (int i = 0; i < M; i++) {
                        sumOfLogLikelihoods += Math.log(n2.density(observations[i]) / n1.density(observations[i]));
                    }
                    double currentScore = Math.max(0.0, previousScore + sumOfLogLikelihoods);
                    previousScore = currentScore;
                    if (currentScore > threshold) {
                        if (iteration < threshForChange) {
                            totalFalseAlarms++;
                        } else if (iteration >= threshold) {
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
            falseAlarmRates.add(ratio);
            detectionDelays.add(averageDetectionDelay);
        }       
        System.out.println("To conclude: false alarm rates and detection delays:");
        System.out.println(falseAlarmRates);
        System.out.println(detectionDelays);
    }
    
    /**
     * This is implemented similar to the single-sensor test, except we have to do things in parallel so we'll
     * have (among other things) a previousScores array that has 10 values.
     * 
     * @param n1
     * @param n2
     * @param M The number of sensors
     */
    public static void testParallelPageOneStateGaussian(NormalDistribution n1, NormalDistribution n2, int M) {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.25) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            double averageDetectionDelay = -1.0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= 10000; trial++) {
                boolean pageTestDone = false;
                int threshForChange = 500;
                int iteration = 0;
                double[] previousScores = new double[M];
                while (!pageTestDone) {
                    double[] observations = new double[M];
                    if (iteration < threshForChange) {
                        totalSamplesFalse++;
                        for (int i = 0; i < M; i++) {
                            observations[i] = n1.sample();
                        }
                    } else {
                        for (int i = 0; i < M; i++) {
                            observations[i] = n2.sample();
                        }
                    }
                    double[] currentScores = new double[M];
                    for (int i = 0; i < M; i++) {
                        double likelihoodComponent = Math.log(n2.density(observations[i]) / n1.density(observations[i]));
                        currentScores[i] = Math.max(0.0, previousScores[i] + likelihoodComponent);
                    }
                    for (int i = 0; i < M; i++) {
                        previousScores[i] = currentScores[i];
                    }
                    // Now we check if any of the M sensors has passed the threshold
                    // Note: since we're only taking the max, we do not "double count" false alarms.
                    double highestScore = maxArrayValue(currentScores);
                    if (highestScore > threshold) {
                        if (iteration < threshForChange) {
                            totalFalseAlarms++;
                        } else if (iteration >= threshold) {
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
            falseAlarmRates.add(ratio);
            detectionDelays.add(averageDetectionDelay);
        }
        System.out.println("To conclude: false alarm rates and detection delays:");
        System.out.println(falseAlarmRates);
        System.out.println(detectionDelays);
    }
    
    /**
     * The purpose of this is to test out a single sensor when the data is generated from a single source,
     * from two different Gaussian distributions, where the threshold that the Gaussian distribution changes
     * is set at the 500^th-indexed iteration throughout all trials. We only end trials after we see a score
     * that exceeds the threshold. False positives are counted when Page's test exceeds the gamma threshold
     * before the 500^th-indexed iteration. Right now this should be done, though we might tweak some values
     * to make the resulting graphs/plots more readable.
     * 
     * @param n1 The benign normal distribution, which is N(0, 1)
     * @param n2 The malignant normal distribution, which is N(0, 1.032)
     */
    public static void testSinglePageOneStateGaussian(NormalDistribution n1, NormalDistribution n2) {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.25) {
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
                    double likelihoodComponent = Math.log(n2.density(observation) / n1.density(observation));
                    double currentScore = Math.max(0.0, previousScore + likelihoodComponent);
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
            falseAlarmRates.add(ratio);
            detectionDelays.add(averageDetectionDelay);
        }
        System.out.println("To conclude: false alarm rates and detection delays:");
        System.out.println(falseAlarmRates);
        System.out.println(detectionDelays);
    }
   
    // Some various utilities

    /** A trivial helper function that takes the maximum element of an array.*/
    private static double maxArrayValue(double[] a) {
        double maxValue = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > maxValue) {
                maxValue = a[i];
            }
        }
        return maxValue;
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
