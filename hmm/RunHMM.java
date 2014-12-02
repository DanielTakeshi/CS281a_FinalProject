package hmm;

import utilities.*;
import java.util.*;

/**
 * This will serve as the starting point of the whole software package. It calls the relevant
 * code to generate a list of observations, which are then analyzed with our HMM and feature
 * tracking machinery so as to determine when the anomalous event occurs (w/high probability).
 * For now, assume we have a fixed observation sequence of some length, and we're going to
 * sample observations from two different "generators" (one of which is the anomalous HMM).
 * 
 * The simple test case flag lets us look at these manually; useful for debugging or small cases.
 * 
 * @author Daniel Seita
 */
public class RunHMM {
    
    //public static final int NUM_OBSERVATIONS = 10;
    //private static List<Observation> observations;

	public static void main(String[] args) {
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
	}
	
	/**
	 * This will go through the states and generate the observation sequence. Ideally, we go through
	 * each state and run our statistical tests on it to determine likelihoods. I'm still kind of fuzzy
	 * on how this will work.
	 * 
	 * @param stateList
	 */
	private static void generateObservations(List<State> stateList) {
	    // TODO Implement
	}
	
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
	private static List<State> mixUpStates(List<State> benignStates, List<State> researchStates) {
	    List<State> result = new ArrayList<State>();
	    // TODO Implement
	    return result;
	}
	
	/**
	 * For now, we'll have this return a list of the hidden states. Then we'll put the functionality
	 * to generate observations in a different method.
	 * 
	 * @param thisHMM A Hidden Markov Model that we can run through.
	 * @return 
	 */
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

	/**
	 * Given a state and a value, draw the next state, which is either itself or from the 
	 * state's children.
	 * 
	 * @param state
	 * @param val
	 * @return
	 */
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

}
