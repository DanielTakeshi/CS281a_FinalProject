package hmm;

import java.util.*;

/**
 * This will serve as the starting point of the whole software package. It calls the relevant code 
 * generate a list of observations, which are then analyzed with our HMM machinery so as to determine
 * when anomalous events start to occur (w/high probability). Start with some simplifying assumptions.
 * 
 * @author Daniel Seita
 */
public class RunHMM {
    
	public static void main(String[] args) {

	    // Some administrative/statistics set up
	    int n = 5;
	    int threshold = 50; // Indicates the splitting point
	    int maxObservations = 150;
	    HMM benign = new HMM(convertToLogProbs(createRandomTransitions(n)),
	            convertToLogProbs(createRandomEmissions(n, false)));
	    HMM adversary = new HMM(convertToLogProbs(createRandomTransitions(n)),
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
	    int currentState = 0;
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
	            currentState = benign.generateNextState(currentState);
	            value = benign.generateObservation(currentState);
	            System.out.println("State = " + currentState + ", observation = " + value);
	        } else {
	            if (obsIndex == threshold) {
	                currentState = 0; // We'll just default to 0 the first time we are in adversary
	            } else {
	                currentState = adversary.generateNextState(currentState);
	            }
	            value = adversary.generateObservation(currentState);
	            System.out.println("State = " + currentState + ", value = " + value);
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
	            if (valueToAdd > 100) {
	                valueToAdd = 100;
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
	            if (valueToAdd > 100) {
	                valueToAdd = 100; // For readability
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
	
	/**
	 * Gets the "clamped" list from the full observations list and the cumulative sums list, so it only
	 * includes the points after the last reset. For example, fullList could be [1,0,1,0,0,1] corresponding
	 * to cumulativeSum = [4.5, 2.1, 0, 0, 3,5, 5.0], so the list to return would be [0,1], i.e., the last
	 * two elements of the full list.
	 * 
	 * @param fullList
	 * @param cumulativeSum
	 * @return
	 */
	public static List<Integer> getClampedList(List<Integer> fullList, List<Double> cumulativeSum) {
	    return null;
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
            Utilities utils = new Utilities();
            double sum = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < a.length; i++) {
                sum = utils.logAdd(sum, a[i]);
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
            Utilities utils = new Utilities();
            double sum = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < a.length; i++) {
                sum = utils.logAdd(sum, a[i]);
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


/**
 * Class for HMMs I'll use. Keep in mind that calculations are to be done in log space, so be sure to convert
 * to normal-space if it makes things easier.
 * 
 * @author Daniel Seita
 */
class HMM {
    
    public int numStates;
    public double[][] transitions;
    public double[][] emissions;
    public Utilities utils;
    
    public HMM(double[][] transitions, double[][] emissions) {
        this.numStates = transitions.length;
        this.transitions = transitions;
        this.emissions = emissions;
        utils = new Utilities();
    }
       
    /**
     * Generates the next state according to the transition matrix. I'll exponentiate to make things easier.
     * 
     * @param state The "previous" state
     * @return The next state, probabilistically chosen
     */
    public int generateNextState(int state) {
        double value = (new Random()).nextDouble();
        double sum = Math.exp(this.transitions[state][0]);
        for (int i = 0; i < numStates-1; i++) {
            if (value < sum) {
                return i;
            }
            sum += Math.exp(this.transitions[state][i+1]);
        }
        return numStates-1;
    }
    
    /**
     * Probabilistically generates a binary-valued observation depending on the state we're in. This is
     * easiest to implement when comparing with the exponentiated version of this.emissions[state][0].
     * 
     * @param state The current state of the HMM
     * @return The observation generated from the state, probabilistically chosen
     */
    public int generateObservation(int state) {
        if ((new Random()).nextDouble() < Math.exp(this.emissions[state][0])) {
            return 0;
        } else {
            return 1;
        }
    }
    
    /**
     * Given a list of observations generated (possibly NOT from this particular HMM), compute the forward
     * probability, i.e., the probability of generating that observation given this HMM. Note: we use log
     * probabilities, and also we only have the list of observations, so we need to sum over all possible
     * current states we could be in for the last observation.
     * 
     * @param observations The list of observations.
     * @return The log of the forward probability.
     */
    public double getForwardProbability(List<Integer> observations) {
        double result = Double.NEGATIVE_INFINITY;
        double[][] tempArray = new double[observations.size()][numStates];
        for (int state = 0; state < numStates; state++) {
            double logProb = this.transitions[0][state] + this.emissions[state][observations.get(0)];
            assert !Double.isNaN(logProb) : "Problem: logProb = " + logProb;
            result = utils.logAdd(result, logProb);
            tempArray[0][state] = logProb;
        }
        assert Math.exp(result) >= 0.0 && Math.exp(result) <= 1.0 : "Problem: Math.exp(result) = " + Math.exp(result);
        if (observations.size() == 1) return result;
        for (int obs = 1; obs < observations.size(); obs++) {
            for (int j = 0; j < numStates; j++) {
                double logSum = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < numStates; i++) {
                    double logComponent = tempArray[obs-1][i] + transitions[i][j] + emissions[j][observations.get(obs)];
                    assert !Double.isNaN(logComponent) : "logComponent = " + logComponent;
                    logSum = utils.logAdd(logSum, logComponent);
                }
                tempArray[obs][j] = logSum;
            }
        }
        result = Double.NEGATIVE_INFINITY;
        for (int state = 0; state < numStates; state++) {
            result = utils.logAdd(result, tempArray[observations.size()-1][state]);
        }
        return result;
    }
}


/**
 * Just some utilities that I'll eventually put in its own file
 * 
 * @author Daniel Seita
 */
class Utilities {
    
    public Utilities() {
        // TODO put something here
    }
    
    /**
     * Given log(x) and log(y), returns log(x+y). This method is from Dan Klein in the SloppyMath
     * class from the Berkeley NLP Parser. This formula is pretty commonplace, though.
     * 
     * @param logX This is log(x) for some x.
     * @param logY This is log(y) for some y.
     * @return The value log(x+y).
     */
	public double logAdd(double logX, double logY) {
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ///////////////////// 
    // OLD STUFF BELOW //
    ///////////////////// 

	    /*
	    // The old main method, let me save this until I really need it
	    System.out.println("*** Now starting up the HMM code ***");
	    boolean simpleTestCase = false;
		Map<String, String> argMap = CommandLineUtils.simpleCommandLineParser(args);
		if (argMap.containsKey("-simple")) {
		    simpleTestCase = true;
		}
		ResearchHiddenMarkovModel researchHMM = new ResearchHiddenMarkovModel();
		BenignHiddenMarkovModel benignHMM = new BenignHiddenMarkovModel(100);
		List<State> benignStates = generateHiddenStateSequence(benignHMM);
		List<State> researchStates = generateHiddenStateSequence(researchHMM);
		List<State> mixedStateSequence = mixUpStates(benignStates, researchStates);
		generateObservations(mixedStateSequence);
		*/
	
	/**
	 * This will go through the states and generate the observation sequence. Ideally, we go through
	 * each state and run our statistical tests on it to determine likelihoods. I'm still kind of fuzzy
	 * on how this will work.
	 * 
	 * @param stateList
	 */
	//private static void generateObservations(List<State> stateList) {
	    // TODO Implement
	//}
	
	/**
	 * Given state sequences from different HMMs, we need to determine how to mix them together to
	 * form one sequence. This will be more of a black art, and needs some randomnes. Ideally, the
	 * benign state sequence is much longer than the other, so we can run it for a while and then
	 * "flip on" the research states so that we can do a, say, 75-25 split? Then after the research
	 * states are done, we should still have some more benign states to run.
	 * 
	 * @param benignStates
	 * @param researchStates
	 * @return
	 */
    /*
	private static List<State> mixUpStates(List<State> benignStates, List<State> researchStates) {
	    List<State> result = new ArrayList<State>();
	    // TODO Implement
	    return result;
	}*/
	
	/**
	 * For now, we'll have this return a list of the hidden states. Then we'll put the functionality
	 * to generate observations in a different method.
	 * 
	 * @param thisHMM A Hidden Markov Model that we can run through.
	 * @return 
	 */
	/*
	private static List<State> generateHiddenStateSequence(HiddenMarkovModel thisHMM) {
	    List<State> stateList = new ArrayList<State>();
	    Random rand = new Random();
	    State state = thisHMM.getRoot();
	    int obsIndex = 0;
	    stateList.add(state);
	    while (!state.isEndState()) {
	        obsIndex++;
	        double val = rand.nextDouble();
	        //System.out.println("\nState number " + obsIndex);
	        //System.out.println("Current state: " + state.toString());
	        //System.out.println("Value drawn was " + val);
	        state = determineNextState(state, val);
	        stateList.add(state);
	    }
	    return stateList;
	}
	*/

	/**
	 * Given a state and a value, draw the next state, which is either itself or from the 
	 * state's children.
	 * 
	 * @param state
	 * @param val
	 * @return
	 */
	/*
	private static State determineNextState(State state, double val) {
	    State nextState = null;
	    boolean newState = false;
	    if (val <= state.getSelfLoopProb()) {
	        nextState = state;
	        newState = true;
	    } else {
	        double[] probThresh = new double[state.numChildren() + 1];
	        probThresh[0] = state.getSelfLoopProb();
	        for (int child = 1; child <= state.numChildren(); child++) {
	            probThresh[child] = probThresh[child-1] + state.getProbabilities().get(child-1);
	        }
	        if (Math.abs(probThresh[probThresh.length-1] - 1.0) > 0.000001) {
	            System.out.println("ERROR: probThresh[n] is not close to one: " + probThresh[probThresh.length-1]);
	            System.exit(0);
	        }
	        for (int i = 1; i < probThresh.length; i++) {
	            if (newState) continue;
	            if (val <= probThresh[i]) {
	                nextState = state.getChildren().get(i-1); // Little confusing
	                newState = true;
	            }
	        }
	    }
	    return nextState;
	}
	*/