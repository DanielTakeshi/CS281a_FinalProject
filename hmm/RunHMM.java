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
 * TODO Extend later to make it really "streaming" instead of getting all the observations at
 * once...
 * 
 * TODO Create an HMM superclass?
 * 
 * @author Daniel Seita
 */
public class RunHMM {
    
    public static final int NUM_OBSERVATIONS = 10;
    //private static List<Observation> observations;

	public static void main(String[] args) {
	    System.out.println("*** Now starting up the HMM code ***");
	    /*
	    boolean verbose = false;
		Map<String, String> argMap = CommandLineUtils.simpleCommandLineParser(args);
		if (argMap.containsKey("-verbose")) {
		    verbose = true;
		}
		*/
		ResearchDevelopmentHMM rdHMM = new ResearchDevelopmentHMM();
		System.out.println("Here's the RD HMM:\n");
		rdHMM.printRD();
		System.out.println();
		generateObservationSequence(rdHMM);
	}
	
	/**
	 * 
	 * @param rdHMM
	 */
	private static void generateObservationSequence(ResearchDevelopmentHMM rdHMM) {
	    Random rand = new Random();
	    State state = rdHMM.root;
	    int obsIndex = 0;
	    while (!state.isEndState()) {
	        obsIndex++;
	        System.out.println("\nObservation Index = " + obsIndex);
	        System.out.println("Current state: " + state.toString());
	        double probability = rand.nextDouble();
	        System.out.println("Value drawn was " + probability);
	        if (probability <= state.getSelfLoopProb()) {
	            continue;
	        }
	        double[] probThresh = new double[state.numChildren() + 1];
	        probThresh[0] = state.getSelfLoopProb();
	        for (int child = 1; child <= state.numChildren(); child++) {
	            probThresh[child] = probThresh[child-1] + state.getProbabilities().get(child-1);
	        }
	        boolean newState = false;
	        for (int i = 1; i < probThresh.length; i++) {
	            if (newState) continue;
	            if (i <= probThresh[i]) {
	                state = state.getChildren().get(i-1); // Little confusing
	                newState = true;
	            }
	        }
	        if(!newState) {
	            System.out.println("ERROR: newState is false");
	            System.exit(0);
	        }
	    }
	    System.out.println("Final state: " + state.toString());
	}

}
