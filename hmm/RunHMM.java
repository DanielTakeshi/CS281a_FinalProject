package hmm;

import edu.berkeley.nlp.util.*;
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
    private static List<Observation> observations;

	public static void main(String[] args) {
		Map<String, String> argMap = CommandLineUtils.simpleCommandLineParser(args);
		// TODO use argMap for parsing command line arguments
		NullHypothesis nullHypothesis = new NullHypothesis();
		ResearchDevelopmentHMM rdHMM = new ResearchDevelopmentHMM();
		observations = generateObservationSequence(nullHypothesis, rdHMM);
		runSequence(observations, nullHypothesis, rdHMM);
	}
	
	/**
	 * 
	 * @param nullHypothesis
	 * @param rdHMM
	 * @return
	 */
	private static List<Observation> generateObservationSequence(NullHypothesis nullHypothesis, ResearchDevelopmentHMM rdHMM) {
	    // TODO implement
	    return null;
	}
	
	/**
	 * 
	 * @param observations
	 * @param nullHypothesis
	 * @param rdHMM
	 */
	private static void runSequence(List<Observation> observations, NullHypothesis nullHypothesis, ResearchDevelopmentHMM rdHMM) {
	    // TODO implement
	}
	
}
