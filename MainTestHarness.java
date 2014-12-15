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
        //System.out.println("Now testing with a one-state, Gaussian source.");
        //conductSingleGaussianSourceTests();
        System.out.println("Done w/one-state Gaussian source. Now doing the HMM Gaussian source.");
        conductHMMGaussianSourceTests();
        System.out.println("Done.");
    }
    
    /** Conducts the HMM Gaussian source tests. Comment out if desired. */
    public static void conductHMMGaussianSourceTests() {
        int numSensors = 10;
        int numHMMStates = 20;
        NormalDistribution benignStdDev = new NormalDistribution(0.0, 1.0); // Used to generate the benign standard deviations
        NormalDistribution[] benignNormals = new NormalDistribution[numHMMStates];
        NormalDistribution malignantStdDev = new NormalDistribution(0.0, 4.0); // Need to make this different from benign
        NormalDistribution[] malignantNormals = new NormalDistribution[numHMMStates];
        // Use a "folded normal" distribution of the standard deviations with Math.abs( ... ).
        for (int i = 0; i < numHMMStates; i++) {
            benignNormals[i] = new NormalDistribution(0, Math.abs(benignStdDev.sample()));
            System.out.println("Benign " + i + " has st dev " + benignNormals[i].getStandardDeviation());
            malignantNormals[i] = new NormalDistribution(0.1, Math.abs(malignantStdDev.sample()));
            System.out.println("Malignant " + i + " has st dev " + malignantNormals[i].getStandardDeviation());
        }
        NormalHiddenMarkovModel hmm1 = new NormalHiddenMarkovModel(numHMMStates, benignNormals);
        NormalHiddenMarkovModel hmm2 = new NormalHiddenMarkovModel(numHMMStates, malignantNormals);
        System.out.println("\n\nNow doing a single Page's test ...");
        testSinglePageHMMGaussian(hmm1, hmm2, numSensors);
        System.out.println("\n\nNow doing Page's test in parallel ...");
        testParallelPageHMMGaussian(hmm1, hmm2, numSensors);
        System.out.println("\n\nNow doing the centralized entity Test ...");
        testIdealCenterHMMGaussian(hmm1, hmm2, numSensors);
        System.out.println("\n\nNow doing running consensus ...");
        testRunningConsensusHMMGaussian(hmm1, hmm2, numSensors);
    }
    
    /** Conducts the single-source, Gaussian tests. Comment out if desired. */
    public static void conductSingleGaussianSourceTests() {
        int numSensors = 10;
        NormalDistribution n1 = new NormalDistribution(0.0, 1.0);
        NormalDistribution n2 = new NormalDistribution(0.0, 1.032); // The doc says 1.032 is stdev, NOT variance
        System.out.println("\n\nNow doing a single Page's test ...");
        testSinglePageOneStateGaussian(n1, n2);
        System.out.println("\n\nNow doing Page's test in parallel ...");
        //testParallelPageOneStateGaussian(n1, n2, numSensors);
        System.out.println("\n\nNow doing the centralized entity Test ...");
        //testIdealCenterOneStateGaussian(n1, n2, numSensors);
        System.out.println("\n\nNow doing running consensus ...");
        //testRunningConsensusOneStateGaussian(n1, n2, numSensors);       
    }
    
    /**
     * Running consensus for the HMM Gaussian case
     * 
     * @param hmm1 The benign HMM
     * @param hmm2 The malignant HMM
     * @param M The number of sensors
     */
    public static void testRunningConsensusHMMGaussian(NormalHiddenMarkovModel hmm1, NormalHiddenMarkovModel hmm2, int M) {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.1) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            double averageDetectionDelay = -1.0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= 10000; trial++) {
                // TODO
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, averageDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add(averageDetectionDelay);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
    }

    /**
     * The ideal center for the HMM Gaussian case
     * 
     * @param hmm1
     * @param hmm2
     * @param M
     */
    public static void testIdealCenterHMMGaussian(NormalHiddenMarkovModel hmm1, NormalHiddenMarkovModel hmm2, int M) {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();       
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.1) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            double averageDetectionDelay = -1.0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= 10000; trial++) {
                // TODO
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, averageDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add(averageDetectionDelay);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
    }

    /**
     * The parallel bank for the HMM Gaussian case
     * 
     * @param hmm1
     * @param hmm2
     * @param M
     */
    public static void testParallelPageHMMGaussian(NormalHiddenMarkovModel hmm1, NormalHiddenMarkovModel hmm2, int M) {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();       
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.1) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            double averageDetectionDelay = -1.0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= 10000; trial++) {
                // TODO
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, averageDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add(averageDetectionDelay);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
    }

    /**
     * The single sensor for the HMM Gaussian case
     * 
     * @param hmm1
     * @param hmm2
     * @param M
     */
    public static void testSinglePageHMMGaussian(NormalHiddenMarkovModel hmm1, NormalHiddenMarkovModel hmm2, int M)  {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();       
        for (double threshold = 1.0; threshold <= 20; threshold += 1.0) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            double averageDetectionDelay = -1.0;
            System.out.println("Currently using threshold gamma = " + threshold);
            
            // For each trial, reset HMM states. Then proceed as normal but need to have list of observations.
            for (int trial = 1; trial <= 5000; trial++) {
                if (trial % 1000 == 0) System.out.println("trial " + trial);
                hmm1.resetState();
                hmm2.resetState();
                boolean pageTestDone = false;
                int threshForChange = 500;
                int iteration = 0;
                List<Double> clampedObservations = new ArrayList<Double>();
                double previousScore = 0.0;
                double[] cachedLogForwardProbabilities1 = new double[hmm1.getNumStates()]; // Cache of log probs
                double[] cachedLogForwardProbabilities2 = new double[hmm2.getNumStates()]; // Cache of log probs
                while (!pageTestDone) {
                    double observation = 0.0;
                    if (iteration < threshForChange) {
                        totalSamplesFalse++;
                        hmm1.generateNextState();
                        observation = hmm1.generateObservation();
                    } else {
                        hmm2.generateNextState();
                        observation = hmm2.generateObservation();
                    }
                    clampedObservations.add(observation);
                    
                    // These take too long
                    //double benignObservation = hmm1.getLogForwardProbability(clampedObservations);
                    //double malignantObservation = hmm2.getLogForwardProbability(clampedObservations);
                    
                    // So let's do it a faster way.
                    if (clampedObservations.size() == 1) {
                        cachedLogForwardProbabilities1 = hmm1.cachedLogForwardProbabilityBaseCase(observation);
                        cachedLogForwardProbabilities2 = hmm2.cachedLogForwardProbabilityBaseCase(observation);
                    } else {
                        cachedLogForwardProbabilities1 = hmm1.cachedLogForwardProbability(
                                cachedLogForwardProbabilities1, clampedObservations);
                        cachedLogForwardProbabilities2 = hmm2.cachedLogForwardProbability(
                                cachedLogForwardProbabilities2, clampedObservations);
                    }
                    double benignObservation = Double.NEGATIVE_INFINITY;
                    for (int i = 0; i < hmm1.getNumStates(); i++) {
                        benignObservation = logSum(benignObservation, cachedLogForwardProbabilities1[i]);
                    }
                    double malignantObservation = Double.NEGATIVE_INFINITY;
                    for (int i = 0; i < hmm2.getNumStates(); i++) {
                        malignantObservation = logSum(malignantObservation, cachedLogForwardProbabilities2[i]);
                    }

                    // Now back to normal...remember that this is in LOG space and that MALIGNANT comes BEFORE the BENIGN
                    double likelihoodComponent = malignantObservation - benignObservation;
                    double currentScore = Math.max(0.0, previousScore + likelihoodComponent);
                    previousScore = currentScore;
                    if (currentScore > threshold) {
                        if (iteration < threshForChange) {
                            totalFalseAlarms++;
                        } else {
                            pageTestDone = true;
                            int thisDetectionDelay = iteration - 500;
                            averageDetectionDelay = ((trial-1)*averageDetectionDelay + thisDetectionDelay) / ((double) trial);
                        }
                    }
                    if (currentScore == 0.0) {
                        clampedObservations = new ArrayList<Double>();  // Reset the observation list
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, averageDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add(averageDetectionDelay);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
    }

    /**
     * This will be perhaps the most interesting part of this code: the Running Consensus scheme. One
     * thing that is not clear from the paper is how they probabilistically set up the v = 5 exchanges.
     * The way I will do it is that I'll just pick five random exchanges and take the average. Note: 
     * some of these may be repeated, e.g., I might choose (1,2) and later get (1,2) again.
     * 
     * @param n1
     * @param n2
     * @param M The number of sensors
     */
    public static void testRunningConsensusOneStateGaussian(NormalDistribution n1, NormalDistribution n2, int M) {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.1) {
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

                // When performing Running Consensus, we will fix the decision to be made at the first sensor
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
                    
                    // Set up the equation [S_n] = W_n*[S_{n-1}] + M*W_n*[l(x)], so first set up the last double[].
                    double[] logLikelihoods = new double[M];
                    for (int i = 0; i < M; i++) { // Include M factor because that's what they did in the paper
                        logLikelihoods[i] = M * Math.log(n2.density(observations[i]) / n1.density(observations[i]));
                    }

                    // Now pick the 5 pairs we need to exchange
                    for (int pair = 0; pair < 5; pair++) {
                        int firstIndex = randInt(0, M-1);
                        int secondIndex = randInt(0, M-1);
                        while (secondIndex == firstIndex) {
                            secondIndex = randInt(0, M-1);
                        }
                        double firstAverage = (previousScores[firstIndex] + previousScores[secondIndex]) / 2.0;
                        previousScores[firstIndex] = firstAverage;
                        previousScores[secondIndex] = firstAverage;
                        double secondAverage = (logLikelihoods[firstIndex] + logLikelihoods[secondIndex]) / 2.0;
                        logLikelihoods[firstIndex] = secondAverage;
                        logLikelihoods[secondIndex] = secondAverage;
                    }
                    
                    // We're done with the pairwise stuff so now it's a matter of adding two things together.
                    double[] currentScores = new double[M];
                    for (int i = 0; i < M; i++) {
                        currentScores[i] = Math.max(0.0, previousScores[i] + logLikelihoods[i]);
                    }
                    for (int i = 0; i < M; i++) {
                        previousScores[i] = currentScores[i];
                    }
                    if (currentScores[0] > threshold) {
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
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, averageDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add(averageDetectionDelay);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
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
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.1) {
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
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, averageDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add(averageDetectionDelay);
        }       
        printFinalStatistics(falseAlarmRates, detectionDelays);
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
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.1) {
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
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, averageDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add(averageDetectionDelay);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
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
        for (double threshold = 1.0; threshold <= 5.0; threshold += 0.1) {
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
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, averageDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add(averageDetectionDelay);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
    }
   
    //////////////////////////////
    // Various Helper Functions //
    //////////////////////////////
    
    /** Prints some statistics after a single threshold's iteration has concluded. Saves some code space. */
    private static void printStatisticsAfterThresh(int totalFalseAlarms, int totalSamplesFalse, int totalSamples,
            double averageDetectionDelay) {
        System.out.println("False alarms = " + totalFalseAlarms);
        System.out.println("Total samples (before thresh) = " + totalSamplesFalse);
        System.out.println("Total samples = " + totalSamples);
        double ratio = ((double) totalFalseAlarms) / totalSamplesFalse;
        System.out.println("Ratio = " + ratio);
        System.out.println("Average detection delay = " + averageDetectionDelay);       
    }

    /** Prints the false alarm & detection delay lists after going through all iterations for a test case */
    private static void printFinalStatistics(List<Double> falseAlarmRates, List<Double> detectionDelays) {
        System.out.println("To conclude: false alarm rates and detection delays:");
        System.out.println(falseAlarmRates);
        System.out.println(detectionDelays);       
    }
    
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
     * Returns a pseudo-random number between min and max, inclusive. The difference between min and
     * max can be at most <code>Integer.MAX_VALUE - 1</code>.
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

    /**
     * Given log(x) and log(y), returns log(x+y). This method is from Dan Klein, but the formula is
     * pretty commonplace.
     * 
     * @param logX This is log(x) for some x.
     * @param logY This is log(y) for some y.
     * @return The value log(x+y).
     */
	public static double logSum(double logX, double logY) {
		if (logY > logX) {
			double temp = logX;
			logX = logY;
			logY = temp;
		}
		if (logX == Double.NEGATIVE_INFINITY) { return logX; }
		double negDiff = logY - logX;
		if (negDiff < -20) { return logX; }
		return logX + java.lang.Math.log(1.0 + java.lang.Math.exp(negDiff));
	}
}
