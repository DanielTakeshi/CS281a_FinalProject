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
    private static final int MAX_PAGE_ITERATIONS = 100000;      // Should not see this (normally)
    private static final double HMM_GAUSSIAN_THRESH = 1.0;      // We're starting with 1.0 and going up to this
    private static final double HMM_GAUSSIAN_INCR = 1.0;        // How much we increment the threshold.
    private static final double HMM_GAUSSIAN_TRIALS = 1000.0;   // For single-source Gaussian this was 10,000
    private static final double SINGLE_GAUSSIAN_TRIALS = 10000.0;   // This is 10,000

    /**
     * This will run a bunch of experiments and print to standard output, so it's probably best to
     * save all this info in random files. When plotting, be sure to use pyplot's y-scale setup if
     * using the false alarm rate, because that depends exponentially on the threshold. May need to
     * fix axes as well by using plt.xlim(start, end) etc.
     * 
     * @param args Nothing for now
     */
    public static void main(String[] args) {
        boolean doSingleSourceGaussian = false;     // Will usually be false from now on
        boolean doHMMSourceGaussian = true;
        boolean doDistributionTest = true;
        if (doSingleSourceGaussian) {
            System.out.println("Now testing with a one-state, Gaussian source.");
            conductSingleGaussianSourceTests();
        }
        if (doHMMSourceGaussian) {
            System.out.println("Now testing with an HMM Gaussian source.");
            conductHMMGaussianSourceTests(doDistributionTest);
        }
        System.out.println("Done.");
    }
    
    /** 
     * Conducts the HMM Gaussian source tests. Comment out if desired. Notice that this will play a key
     * role in how well our results look like. It's kind of a black art to figure out the best set of
     * HMM states since it's so easy to get no false alarms or all false alarms, which doesn't give us
     * info on detection delay.
     * 
     * @param doDistributionTest If true, we do a special distribution test designed to check if we are
     *      going to do all four tests or if we're only doing the single page's test multiple times.
     */
    public static void conductHMMGaussianSourceTests(boolean doDistributionTest) {
        int numSensors = 10;
        int numHMMStates = 25;
        NormalDistribution benignStdDev = new NormalDistribution(0.0, 1.0);     // Generates benign s.d.
        NormalDistribution malignantStdDev = new NormalDistribution(0.0, 4.35); // Seems like 4-ish is good?
        System.out.println("Here, detection tests are from HMMs with " + numHMMStates + " states.");
        if (doDistributionTest) {
            List<Double> rates = new ArrayList<Double>();
            List<Double> ratesWithoutMaxIterations = new ArrayList<Double>();
            for (int k = 0; k < 100; k++) {
                System.out.println("On test " + k);
                NormalDistribution[] benignNormals = new NormalDistribution[numHMMStates];
                NormalDistribution[] malignantNormals = new NormalDistribution[numHMMStates];
                for (int i = 0; i < numHMMStates; i++) {
                    benignNormals[i] = new NormalDistribution(0, Math.abs(benignStdDev.sample())); // Folded normal
                    malignantNormals[i] = new NormalDistribution(0.1, Math.abs(malignantStdDev.sample()));
                }
                NormalHiddenMarkovModel hmm1 = new NormalHiddenMarkovModel(numHMMStates, benignNormals);
                NormalHiddenMarkovModel hmm2 = new NormalHiddenMarkovModel(numHMMStates, malignantNormals);
                double falseAlarmRate = testSinglePageHMMGaussian(hmm1, hmm2, numSensors);               
                rates.add(falseAlarmRate); 
                if (falseAlarmRate >= 0.0) {
                    ratesWithoutMaxIterations.add(falseAlarmRate);
                }
                // Note: I changed the method to return -1 if maxIterations has occurred (distributions are too similar)
            }
            System.out.println("Done with single Page's test. Here are the f.a. rates:");
            System.out.println(rates);
            System.out.println("And the ones WITHOUT -1, i.e., max iterations. Length = " + ratesWithoutMaxIterations.size());
            System.out.println(ratesWithoutMaxIterations);
        } else {
            NormalDistribution[] benignNormals = new NormalDistribution[numHMMStates];
            NormalDistribution[] malignantNormals = new NormalDistribution[numHMMStates];
            for (int i = 0; i < numHMMStates; i++) {
                benignNormals[i] = new NormalDistribution(0, Math.abs(benignStdDev.sample())); // Folded normal
                malignantNormals[i] = new NormalDistribution(0.1, Math.abs(malignantStdDev.sample()));
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
    }
    
    /** Conducts the single-source, Gaussian tests. Comment out if desired. */
    public static void conductSingleGaussianSourceTests() {
        int numSensors = 10;
        NormalDistribution n1 = new NormalDistribution(0.0, 1.0);
        NormalDistribution n2 = new NormalDistribution(0.0, 1.032); // The doc says 1.032 is stdev, NOT variance
        System.out.println("\n\nNow doing a single Page's test ...");
        testSinglePageOneStateGaussian(n1, n2);
        System.out.println("\n\nNow doing Page's test in parallel ...");
        testParallelPageOneStateGaussian(n1, n2, numSensors);
        System.out.println("\n\nNow doing the centralized entity Test ...");
        testIdealCenterOneStateGaussian(n1, n2, numSensors);
        System.out.println("\n\nNow doing running consensus ...");
        testRunningConsensusOneStateGaussian(n1, n2, numSensors);       
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
        for (double threshold = 1.0; threshold <= HMM_GAUSSIAN_THRESH; threshold += HMM_GAUSSIAN_INCR) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            int totalDetectionDelay = 0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= HMM_GAUSSIAN_TRIALS; trial++) {
                if (trial % 250 == 0) System.out.println("trial " + trial);
                hmm1.resetState();
                hmm2.resetState();
                boolean pageTestDone = false;
                int threshForChange = 500;
                int iteration = 0;

                 // Each of the M elements here has its own "clamped observation" list. Also need to set up cache.
                ArrayList<ArrayList<Double>> listOfClampedObservations = new ArrayList<ArrayList<Double>>();
                ArrayList<double[]> cachedLogForwardProbs1 = new ArrayList<double[]>();
                ArrayList<double[]> cachedLogForwardProbs2 = new ArrayList<double[]>();
                for (int i = 0; i < M; i++) {
                    ArrayList<Double> emptyList = new ArrayList<Double>();
                    listOfClampedObservations.add(emptyList);
                    cachedLogForwardProbs1.add(new double[hmm1.getNumStates()]);
                    cachedLogForwardProbs2.add(new double[hmm2.getNumStates()]);
                }               
                double[] previousScores = new double[M];
                while (!pageTestDone) {
                    double[] observations = new double[M];
                    if (iteration < threshForChange) {
                        totalSamplesFalse++;
                        hmm1.generateNextState();
                        for (int i = 0; i < M; i++) {
                            observations[i] = hmm1.generateObservation();
                        }
                    } else {
                        hmm2.generateNextState();
                        for (int i = 0; i < M; i++) {
                            observations[i] = hmm2.generateObservation();
                        }
                    }               

                    // Set up the equation [S_n] = W_n*[S_{n-1}] + M*W_n*[l(x)], so first set up the last double[].
                    double[] logLikelihoods = new double[M];
                    for (int i = 0; i < M; i++) {
                        listOfClampedObservations.get(i).add(observations[i]);  // Don't forget to add observations!
                        if (listOfClampedObservations.get(i).size() == 1) {
                            cachedLogForwardProbs1.set(i, hmm1.cachedLogForwardProbabilityBaseCase(observations[i]));
                            cachedLogForwardProbs2.set(i, hmm2.cachedLogForwardProbabilityBaseCase(observations[i]));
                        } else {
                            cachedLogForwardProbs1.set(i, hmm1.cachedLogForwardProbability(cachedLogForwardProbs1.get(i), observations[i]));
                            cachedLogForwardProbs2.set(i, hmm2.cachedLogForwardProbability(cachedLogForwardProbs2.get(i), observations[i]));
                        }
                        double benignObservation = Double.NEGATIVE_INFINITY;
                        for (int state = 0; state < hmm1.getNumStates(); state++) {
                            benignObservation = logSum(benignObservation, cachedLogForwardProbs1.get(i)[state]);
                        }
                        double malignantObservation = Double.NEGATIVE_INFINITY;
                        for (int state = 0; state < hmm2.getNumStates(); state++) {
                            malignantObservation = logSum(malignantObservation, cachedLogForwardProbs2.get(i)[state]);
                        }
                        //logLikelihoods[i] = M * Math.log(n2.density(observations[i]) / n1.density(observations[i]));
                        logLikelihoods[i] = M * (malignantObservation - benignObservation); // I really hope this is correct.
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
                    if (currentScores[0] > threshold || iteration >= MAX_PAGE_ITERATIONS) {
                        if (iteration < threshForChange) {
                            totalFalseAlarms++;
                        } else if (iteration >= threshForChange) {
                            if (iteration >= MAX_PAGE_ITERATIONS) System.out.println("At maximum number of iterations.");
                            pageTestDone = true;
                            totalDetectionDelay += (iteration = 500);
                        }
                    }
                    
                    // I think this resetting will work
                    for (int i = 0; i < M; i++) {
                        if (currentScores[i] == 0.0) {
                            listOfClampedObservations.set(i, new ArrayList<Double>());  // Reset the observation list
                        }
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, totalDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add( ((double) totalDetectionDelay) / HMM_GAUSSIAN_TRIALS);
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
        for (double threshold = 1.0; threshold <= HMM_GAUSSIAN_THRESH; threshold += HMM_GAUSSIAN_INCR) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            int totalDetectionDelay = 0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= HMM_GAUSSIAN_TRIALS; trial++) {
                if (trial % 250 == 0) System.out.println("trial " + trial);
                hmm1.resetState();
                hmm2.resetState();
                boolean pageTestDone = false;
                int threshForChange = 500;
                int iteration = 0;

                // Each of the M elements here has its own "clamped observation" list. Also need to set up cache.
                ArrayList<ArrayList<Double>> listOfClampedObservations = new ArrayList<ArrayList<Double>>();
                ArrayList<double[]> cachedLogForwardProbs1 = new ArrayList<double[]>();
                ArrayList<double[]> cachedLogForwardProbs2 = new ArrayList<double[]>();
                for (int i = 0; i < M; i++) {
                    ArrayList<Double> emptyList = new ArrayList<Double>();
                    listOfClampedObservations.add(emptyList);
                    cachedLogForwardProbs1.add(new double[hmm1.getNumStates()]);
                    cachedLogForwardProbs2.add(new double[hmm2.getNumStates()]);
                }
                double previousScore = 0.0; // Must use a single value for the score portion
                
                // Now proceed with Page's test iteration
                while (!pageTestDone) {
                    double[] observations = new double[M];
                    if (iteration < threshForChange) {
                        totalSamplesFalse++;
                        hmm1.generateNextState();
                        for (int i = 0; i < M; i++) {
                            observations[i] = hmm1.generateObservation();
                        }
                    } else {
                        hmm2.generateNextState();
                        for (int i = 0; i < M; i++) {
                            observations[i] = hmm2.generateObservation();
                        }
                    }
                    double sumOfLogRatios = 0.0;
                    for (int i = 0; i < M; i++) {
                        listOfClampedObservations.get(i).add(observations[i]);
                        if (listOfClampedObservations.get(i).size() == 1) {
                            cachedLogForwardProbs1.set(i, hmm1.cachedLogForwardProbabilityBaseCase(observations[i]));
                            cachedLogForwardProbs2.set(i, hmm2.cachedLogForwardProbabilityBaseCase(observations[i]));
                        } else {
                            cachedLogForwardProbs1.set(i, hmm1.cachedLogForwardProbability(cachedLogForwardProbs1.get(i), observations[i]));
                            cachedLogForwardProbs2.set(i, hmm2.cachedLogForwardProbability(cachedLogForwardProbs2.get(i), observations[i]));
                        }
                        double benignObservation = Double.NEGATIVE_INFINITY;
                        for (int state = 0; state < hmm1.getNumStates(); state++) {
                            benignObservation = logSum(benignObservation, cachedLogForwardProbs1.get(i)[state]);
                        }
                        double malignantObservation = Double.NEGATIVE_INFINITY;
                        for (int state = 0; state < hmm2.getNumStates(); state++) {
                            malignantObservation = logSum(malignantObservation, cachedLogForwardProbs2.get(i)[state]);
                        }
                        double likelihoodComponent = malignantObservation - benignObservation;
                        sumOfLogRatios += likelihoodComponent;
                    }

                    double currentScore = Math.max(0.0, previousScore + sumOfLogRatios);
                    previousScore = currentScore;
                    if (currentScore > threshold || iteration >= MAX_PAGE_ITERATIONS) {
                        if (iteration < threshForChange) {
                            totalFalseAlarms++;
                        } else {
                            if (iteration >= MAX_PAGE_ITERATIONS) System.out.println("At maximum number of iterations.");
                            pageTestDone = true;
                            totalDetectionDelay += (iteration - 500);
                        }
                    }
                    if (currentScore == 0) {
                        for (int i = 0; i < M; i++) {
                            listOfClampedObservations.set(i, new ArrayList<Double>());  // Reset the observation list
                        }
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, totalDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add( ((double) totalDetectionDelay) / HMM_GAUSSIAN_TRIALS);
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
        for (double threshold = 1.0; threshold <= HMM_GAUSSIAN_THRESH; threshold += HMM_GAUSSIAN_INCR) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            int totalDetectionDelay = 0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= HMM_GAUSSIAN_TRIALS; trial++) {
                if (trial % 250 == 0) System.out.println("trial " + trial);
                hmm1.resetState();
                hmm2.resetState();
                boolean pageTestDone = false;
                int threshForChange = 500;
                int iteration = 0;

                // Each of the M elements here has its own "clamped observation" list. Also need to set up cache.
                ArrayList<ArrayList<Double>> listOfClampedObservations = new ArrayList<ArrayList<Double>>();
                ArrayList<double[]> cachedLogForwardProbs1 = new ArrayList<double[]>();
                ArrayList<double[]> cachedLogForwardProbs2 = new ArrayList<double[]>();
                for (int i = 0; i < M; i++) {
                    ArrayList<Double> emptyList = new ArrayList<Double>();
                    listOfClampedObservations.add(emptyList);
                    cachedLogForwardProbs1.add(new double[hmm1.getNumStates()]);
                    cachedLogForwardProbs2.add(new double[hmm2.getNumStates()]);
                }
                double[] previousScores = new double[M];
                
                // Now proceed with Page's test iteration
                while (!pageTestDone) {
                    double[] observation = new double[M];
                    if (iteration < threshForChange) {
                        totalSamplesFalse++;
                        hmm1.generateNextState();
                        for (int i = 0; i < M; i++) {
                            observation[i] = hmm1.generateObservation();
                        }
                    } else {
                        hmm2.generateNextState();
                        for (int i = 0; i < M; i++) {
                            observation[i] = hmm2.generateObservation();
                        }
                    }
                    double[] currentScores = new double[M];
                    for (int i = 0; i < M; i++) {
                        listOfClampedObservations.get(i).add(observation[i]);
                        if (listOfClampedObservations.get(i).size() == 1) {
                            cachedLogForwardProbs1.set(i, hmm1.cachedLogForwardProbabilityBaseCase(observation[i]));
                            cachedLogForwardProbs2.set(i, hmm2.cachedLogForwardProbabilityBaseCase(observation[i]));
                        } else {
                            // Note: listOfClampedObservations.get(i) is like List<Double> observations, so only need last one.
                            cachedLogForwardProbs1.set(i, hmm1.cachedLogForwardProbability(cachedLogForwardProbs1.get(i), observation[i]));
                            cachedLogForwardProbs2.set(i, hmm2.cachedLogForwardProbability(cachedLogForwardProbs2.get(i), observation[i]));
                        }
                        double benignObservation = Double.NEGATIVE_INFINITY;
                        for (int state = 0; state < hmm1.getNumStates(); state++) {
                            benignObservation = logSum(benignObservation, cachedLogForwardProbs1.get(i)[state]);
                        }
                        double malignantObservation = Double.NEGATIVE_INFINITY;
                        for (int state = 0; state < hmm2.getNumStates(); state++) {
                            malignantObservation = logSum(malignantObservation, cachedLogForwardProbs2.get(i)[state]);
                        }
                        double likelihoodComponent = malignantObservation - benignObservation;
                        double currentScore = Math.max(0.0, previousScores[i] + likelihoodComponent);
                        currentScores[i] = currentScore;
                    }

                    // Now with the scores, update previous score and check if it exceeds our threshold.
                    for (int i = 0; i < M; i++) {
                        previousScores[i] = currentScores[i];
                    }
                    double highestCurrentScore = maxArrayValue(currentScores);
                    if (highestCurrentScore > threshold || iteration >= MAX_PAGE_ITERATIONS) {
                        if (iteration < threshForChange) {
                            totalFalseAlarms++;
                        } else {
                            if (iteration >= MAX_PAGE_ITERATIONS) System.out.println("Note: we're at max iterations");
                            pageTestDone = true;
                            totalDetectionDelay += (iteration - 500);
                        }
                    }
                    
                    // Now need to potentially reset any observation sequences
                    for (int i = 0; i < M; i++) {
                        if (currentScores[i] == 0.0) {
                            listOfClampedObservations.set(i, new ArrayList<Double>());  // Reset the observation list
                        }
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, totalDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add( ((double) totalDetectionDelay) / HMM_GAUSSIAN_TRIALS);
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
    public static double testSinglePageHMMGaussian(NormalHiddenMarkovModel hmm1, NormalHiddenMarkovModel hmm2, int M)  {
        List<Double> falseAlarmRates = new ArrayList<Double>();
        List<Double> detectionDelays = new ArrayList<Double>();       
        for (double threshold = 1.0; threshold <= HMM_GAUSSIAN_THRESH; threshold += HMM_GAUSSIAN_INCR) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            int totalDetectionDelay = 0;
            System.out.println("Currently using threshold gamma = " + threshold);
            
            // For each trial, reset HMM states. Then proceed as normal but need to have list of observations.
            for (int trial = 1; trial <= HMM_GAUSSIAN_TRIALS; trial++) {
                if (trial % 250 == 0) System.out.println("trial " + trial);
                hmm1.resetState();
                hmm2.resetState();
                boolean pageTestDone = false;
                int threshForChange = 500;
                int iteration = 0;
                List<Double> clampedObservations = new ArrayList<Double>(); // List of observations AFTER last reset
                double previousScore = 0.0;
                double[] cachedLogForwardProbabilities1 = new double[hmm1.getNumStates()]; // Cache of log probs
                double[] cachedLogForwardProbabilities2 = new double[hmm2.getNumStates()]; // Cache of log probs
                while (!pageTestDone) {
                    double observation = 0.0;
                    if (iteration < threshForChange) {
                        totalSamplesFalse++;
                        hmm1.generateNextState();
                        observation = hmm1.generateObservation();
                        //System.out.println("Observation from hmm1 " + observation);
                    } else {
                        hmm2.generateNextState();
                        observation = hmm2.generateObservation();
                        //System.out.println("Observation from hmm2 " + observation);
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
                        cachedLogForwardProbabilities1 = hmm1.cachedLogForwardProbability(cachedLogForwardProbabilities1, observation);
                        cachedLogForwardProbabilities2 = hmm2.cachedLogForwardProbability(cachedLogForwardProbabilities2, observation);
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
                    //System.out.println("iteration " + iteration + " with score " + currentScore);
                    previousScore = currentScore;
                    if (currentScore > threshold || iteration >= MAX_PAGE_ITERATIONS) {
                        if (iteration < threshForChange) {
                            totalFalseAlarms++;
                        } else {
                            if (iteration >= MAX_PAGE_ITERATIONS) {
                                System.out.println("Note: we're at max iterations. Let's break out of this and "
                                        + "return -1 for false alarms because these distributions are too similar.");
                                return -1.0;
                            }
                            pageTestDone = true;
                            totalDetectionDelay += (iteration - 500);
                            //System.out.println("total detection delay = " + totalDetectionDelay);
                        }
                    }
                    if (currentScore == 0.0) {
                        clampedObservations = new ArrayList<Double>();  // Reset the observation list
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, totalDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add( ((double) totalDetectionDelay) / HMM_GAUSSIAN_TRIALS);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
        // New
        return falseAlarmRates.get(0);
    }
    

    ////////////////////////////
    ////////////////////////////
    // Single-source Gaussian //
    ////////////////////////////
    ////////////////////////////


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
        for (double threshold = 0.5; threshold <= 5.0; threshold += 0.1) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            int totalDetectionDelay = 0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= SINGLE_GAUSSIAN_TRIALS; trial++) {
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
                            totalDetectionDelay += (iteration - 500);
                        }
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, totalDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add( ((double) totalDetectionDelay) / SINGLE_GAUSSIAN_TRIALS);
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
        for (double threshold = 0.5; threshold <= 5.0; threshold += 0.1) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            int totalDetectionDelay = 0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= SINGLE_GAUSSIAN_TRIALS; trial++) {
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
                            totalDetectionDelay += (iteration - 500);
                        }
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, totalDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add( ((double) totalDetectionDelay) / SINGLE_GAUSSIAN_TRIALS);
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
        for (double threshold = 0.5; threshold <= 5.0; threshold += 0.1) {
            int totalFalseAlarms = 0;
            int totalSamplesFalse = 0;
            int totalSamples = 0;
            int totalDetectionDelay = 0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= SINGLE_GAUSSIAN_TRIALS; trial++) {
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
                            totalDetectionDelay += (iteration - 500);
                        }
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, totalDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add( ((double) totalDetectionDelay) / SINGLE_GAUSSIAN_TRIALS);
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
        for (double threshold = 0.5; threshold <= 5.0; threshold += 0.1) {
            int totalFalseAlarms = 0;   // Total times Page's test exceeds thresh when it's benign
            int totalSamplesFalse = 0;  // Total samples drawn before malignant distribution begins
            int totalSamples = 0;       // Total samples drawn overall across all trials here
            int totalDetectionDelay = 0;
            System.out.println("Currently using threshold gamma = " + threshold);
            for (int trial = 1; trial <= SINGLE_GAUSSIAN_TRIALS; trial++) {
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
                            totalDetectionDelay += (iteration - 500);
                        }
                    }
                    totalSamples++;
                    iteration++;
                }
            }
            printStatisticsAfterThresh(totalFalseAlarms, totalSamplesFalse, totalSamples, totalDetectionDelay);
            falseAlarmRates.add( ((double) totalFalseAlarms) / totalSamplesFalse );
            detectionDelays.add( ((double) totalDetectionDelay) / SINGLE_GAUSSIAN_TRIALS);
        }
        printFinalStatistics(falseAlarmRates, detectionDelays);
    }
   
    //////////////////////////////
    // Various Helper Functions //
    //////////////////////////////
    
    /** Prints some statistics after a single threshold's iteration has concluded. Saves some code space. */
    private static void printStatisticsAfterThresh(int totalFalseAlarms, int totalSamplesFalse, int totalSamples,
            int totalDetectionDelay) {
        System.out.println("False alarms = " + totalFalseAlarms);
        System.out.println("Total samples (before thresh) = " + totalSamplesFalse);
        System.out.println("Total samples = " + totalSamples);
        double ratio = ((double) totalFalseAlarms) / totalSamplesFalse;
        System.out.println("Ratio = " + ratio);
        System.out.println("Total detection delay = " + totalDetectionDelay);       
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
