package hmm;

import java.util.*;
import utilities.*;

/**
 * This will serve as the starting point of the whole software package. It calls the relevant code to
 * generate a list of observations, which are then analyzed with our HMM machinery so as to determine
 * when anomalous events start to occur (w/high probability). Start with some simplifying assumptions.
 * 
 * @author Daniel Seita
 */
public class MainHarness {
    
    /**
     * Runs a test suite right now. Originally, I had two HMMs, and we'd go through the first (benign)
     * one, then after it ended, we'd start at state 0 again in the adversary HMM and run through the
     * adversary. But how do we implement a mixture of HMMs? What if each HMM stored its internal state
     * and each time we called "generate next state," it would generate the next HMM state based on its
     * internal state? That seems like a cleaner way of doing this. Now yes we will need to run tests
     * while only using the observations so it's kind of clunky because we have to take probabilities
     * of full observation sequences, but maybe that's an interesing question.
     * 
     * @param args
     */
	public static void main(String[] args) {

	    // Some administrative/statistics set up
	    int n = 5;
	    double hThreshold = Double.POSITIVE_INFINITY; // Set this to be something if I want to limit things
	    int threshold = 100; // Indicates the splitting point (UPDATE this time it will be for mixtures)
	    int maxObservations = 500;
	    HiddenMarkovModel benign = new HiddenMarkovModel(convertToLogProbs(createRandomTransitions(n)),
	            convertToLogProbs(createRandomEmissions(n, false)));
	    HiddenMarkovModel adversary = new HiddenMarkovModel(convertToLogProbs(createRandomTransitions(n)),
	            convertToLogProbs(createRandomEmissions(n, true)));
	    System.out.println("Running this example with " + n + " states.");
	    System.out.println("\nHere are the benign tranisitions and emissions:");
	    printMatrix(benign.transitions);
	    System.out.println();
	    printMatrix(benign.emissions);
	    System.out.println("\nAnd here are the adversary transitions and emissions:");
	    printMatrix(adversary.transitions);
	    System.out.println();
	    printMatrix(adversary.emissions);
	    System.out.println("\n\nNow let's go through the HMM.");
	    int value = -1;
	    int numZeros = 0;
	    int numOnes = 0;
	    List<Double> cumulativeSum = new ArrayList<Double>();
	    List<Double> clampedCumulativeSum = new ArrayList<Double>();
	    List<Integer> clampedObservations = new ArrayList<Integer>();
	    int firstIndexToUse = 0;
	    
	    // Now go through the HMM
	    List<Integer> observationList = new ArrayList<Integer>();
	    for (int obsIndex = 0; obsIndex < maxObservations; obsIndex++) {
	        if (obsIndex == firstIndexToUse) {
	            clampedObservations = new ArrayList<Integer>();
	        }
	        System.out.println("\nObservation Index = " + obsIndex);
	        if (obsIndex < threshold) {
	            benign.generateNextState();
	            value = benign.generateObservation();
	            System.out.println("State = " + benign.getCurrentState() + ", observation = " + value);
	        } else {
	            // Old way was that I just used the adversary state. New way, I think I should use a mix.
	            Random rand = new Random();
	            if (rand.nextDouble() < 0.5) {
	                adversary.generateNextState();
	                value = adversary.generateObservation();
	                System.out.println("State = " + adversary.getCurrentState() + ", value = " + value);
	            } else {
	                benign.generateNextState();
	                value = benign.generateObservation();
	                System.out.println("State = " + benign.getCurrentState() + ", value = " + value);
	            }
	        }
	        if (value == 0) numZeros++;
	        if (value == 1) numOnes++;

	        // The above info shouldn't change based on whether we're using clamped vs full lists. First do full:
	        observationList.add(value);
	        System.out.println("Current (0s, 1s) tally = (" + numZeros + ", " + numOnes + ").");
	        double benignForwardLogProb = benign.getForwardProbability(observationList);
	        double adversaryForwardLogProb = adversary.getForwardProbability(observationList);
	        System.out.println("Benign forward log prob: " + benignForwardLogProb);
	        System.out.println("Adversary forward log prob: " + adversaryForwardLogProb);
	        double thisCumulativeSum = adversaryForwardLogProb - benignForwardLogProb;
	        if (cumulativeSum.size() == 0) {
	            cumulativeSum.add(Math.max(0.0, thisCumulativeSum));
	        } else {
	            double valueToAdd = Math.max(0.0, thisCumulativeSum + cumulativeSum.get(cumulativeSum.size()-1));
	            if (valueToAdd > hThreshold) {
	                valueToAdd = hThreshold;
	            }
	            cumulativeSum.add(valueToAdd);
	        }
	        System.out.println("Cumulative observations: " + observationList);
	        System.out.println(cumulativeSum);

	        // Now do clamped version
	        clampedObservations.add(value);
	        benignForwardLogProb = benign.getForwardProbability(clampedObservations);
	        adversaryForwardLogProb = adversary.getForwardProbability(clampedObservations);
	        System.out.println("Benign CLAMPED forward log prob: " + benignForwardLogProb);
	        System.out.println("Adversary CLAMPED forward log prob: " + adversaryForwardLogProb);
	        thisCumulativeSum = adversaryForwardLogProb - benignForwardLogProb;
	        if (clampedCumulativeSum.size() == 0) {
	            clampedCumulativeSum.add(Math.max(0.0, thisCumulativeSum));
	        } else {
	            double valueToAdd = Math.max(0.0, thisCumulativeSum + clampedCumulativeSum.get(clampedCumulativeSum.size()-1));
	            if (valueToAdd > hThreshold) {
	                valueToAdd = hThreshold; // For readability in graphs (so values don't get too large); can also serve as our h threshold.
	            }
	            clampedCumulativeSum.add(valueToAdd); 
	        }
	        if (clampedCumulativeSum.get(clampedCumulativeSum.size()-1) == 0.0) {
	            firstIndexToUse = obsIndex+1;
	        }
	        System.out.println("Clamped observations: " + clampedObservations);
	        System.out.println(clampedCumulativeSum);
	    }
	}
	
	public static void runTestSuite() {
	    // TODO Use this method to run a test suite
	}
	
	/**
	 * Given the number of states, creates a transition matrix 'transitions' where transitions[k] is the
	 * vector of transition probabilities of the form a_{kx} where x is the index inside transitions[k],
	 * so k is the prior state and x is the current one. Right now it initializes values to be uniformly
	 * distributed between 0 and 1, and then normalizes by row. 
	 * 
	 * @param numStates The number of states in the HMM, i.e., transitions.length and transitions[k].length.
	 * @return A random transition matrix, which could be used as the "A" matrix for a null/benign HMM.
	 */
	public static double[][] createRandomTransitions(int numStates) {
	    double[][] transitions = new double[numStates][numStates];
	    Random rand = new Random();
	    for (int row = 0; row < numStates; row++) {
	        double[] nonNormalizedValues = new double[numStates];
	        double sum = 0.0;
	        for (int i = 0; i < numStates; i++) {
	            nonNormalizedValues[i] = rand.nextDouble();
	            sum += nonNormalizedValues[i];
	        }
	        for (int i = 0; i < numStates; i++) {
	            transitions[row][i] = nonNormalizedValues[i] / sum;
	            assert transitions[row][i] >= 0.0 && transitions[row][i] <= 1.0 : "Problem:"
	                    + " transitions[row][i] = " + transitions[row][i];
	        }
	    }
	    return transitions;
	}
	
	/**
	 * Given the number of states, creates the emissions array "emissions" for a given HMM. Here, we'll
	 * assume we have outputs 0 and 1, and that emissions[0] = [Pr(0 | S0), Pr(1 | S0)], so we're going
	 * to condition emissions[k] on the STATE number. Values are initialized randomly. TODO I'm going to
	 * add a special change to if we have adversary vs nonadversary to make things a little more interesting.
	 * Adversary will have a higher chance of generating ones.
	 * 
	 * @param numStates The number of states in the HMM.
	 * @param adversary True if it's the adversary, false if otherwise.
	 * @return A randomly initialized emissions probability matrix (or array of arrays...).
	 */
	public static double[][] createRandomEmissions(int numStates, boolean adversary) {
	    double[][] emissions = new double[numStates][2];
	    Random rand = new Random();
	    for (int state = 0; state < numStates; state++) {
	        double value = rand.nextDouble();
	        if ((value < 0.5 && adversary) || (value > 0.5 && !adversary)) {
	            emissions[state][0] = value;
	            emissions[state][1] = 1.0 - value;
	        } else {
		        emissions[state][1] = value;
	            emissions[state][0] = 1.0 - value;
	        }
	    }
	    return emissions;
	}

    /**
     * A simple method to normalize an existing matrix by ROWS, i.e., as if we were re-normalizing the
     * transition matrix of an HMM after an adversary targets it. It can also be applied to emissions.
     * It assumes that existingMatrix[k].length is the same for all k \in { 0, 1, ... , numStates-1 }.
     * 
     * @param existingMatrix The matrix we need to normalize (it's fine if it's already normalized)
     * @param logSpace True if we're doing this in log space, false if not
     * @return A normalized version of existingMatrix.
     */
    public static double[][] normalizeExistingMatrix(double[][] existingMatrix, boolean logSpace) {
        int numStates = existingMatrix.length;
        double[][] updatedMatrix = new double[existingMatrix.length][existingMatrix[0].length];
        for (int row = 0; row < numStates; row++) {
            if (!checkIfArrayNormalizes(existingMatrix[row], logSpace)) {
                updatedMatrix[row] = normalizeExistingArray(existingMatrix[row], logSpace);
            }
        }
        return updatedMatrix;
    }
    
    /**
     * A simple method to check if an array of doubles is normalized.
     * 
     * @param a The array of doubles that we wish to be normalized.
     * @param logSpace True if we're doing this in log space, false if not
     * @return Whether | sum-1.0 | is less than a tolerance (true) or not (false, doesn't normalize).
     */
    public static boolean checkIfArrayNormalizes(double[] a, boolean logSpace) {
        if (logSpace) {
            double sum = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < a.length; i++) {
                sum = DanielUtilities.logAdd(sum, a[i]);
            }
            return Math.abs(Math.exp(sum) - 1.0) < 0.000001;
        } else {
            double sum = 0.0;
            for (int i = 0; i < a.length; i++) {
                sum += a[i];
            }
            return Math.abs(sum - 1.0) < 0.000001;
        }
    }

    /**
     * A simple method to re-normalize an existing array of doubles.
     * 
     * @param a The array we need to re-normalize (it's fine if it's already normalized, though)
     * @param logSpace True if we're doing this in log space, false if not
     * @return The normalized version of a.
     */
    public static double[] normalizeExistingArray(double[] a, boolean logSpace) {
        double[] updatedArray = new double[a.length];
        if (logSpace) {
            double sum = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < a.length; i++) {
                sum = DanielUtilities.logAdd(sum, a[i]);
            }
            for (int i = 0; i < a.length; i++) {
                updatedArray[i] = a[i] - sum;
                assert Math.exp(updatedArray[i]) >= 0.0 && Math.exp(updatedArray[i]) <= 1.0 : "Problem: "
                        + "Math.exp(updatedArray[i]) = " + Math.exp(updatedArray[i]);
            }
        } else {
            double sum = 0.0;
            for (int i = 0; i < a.length; i++) {
                sum += a[i];
            }
            for (int i = 0; i < a.length; i++) {
                updatedArray[i] = a[i] / sum;   
                assert updatedArray[i] >= 0.0 && updatedArray[i] <= 1.0 : "Problem: updatedArray[i] = " + updatedArray[i];
            }
        }
        return updatedArray;
    }

    /**
     * Useful to print out matrices of doubles, e.g., transition matrix probabilities, and even emissions.
     * 
     * @param matrix The matrix (or if you prefer it this way, an array of arrays) that we want to print.
     */
    public static void printMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            System.out.println(Arrays.toString(matrix[i]));
        }
    }	

    /**
     * Converts a matrix of probabilities to a matrix of log probabilities. This should work as long as the
     * input has the same a[k].length for all k.
     * 
     * @param a 
     * @return
     */
    public static double[][] convertToLogProbs(double[][] a) {
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                result[i][j] = Math.log(a[i][j]);
            }
        }
        return result;
    }
 
}
