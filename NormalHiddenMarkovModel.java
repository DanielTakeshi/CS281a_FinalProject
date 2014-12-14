import java.util.*;
import org.apache.commons.math3.distribution.*;

/**
 * Class for HMMs I'll use. Keep in mind that calculations are to be done in log space, so be sure
 * to convert to normal-space if it makes things easier. EDIT: this class assumes that the emissions
 * are Gaussian. To make things different (e.g., Bernoulli), use a different class.
 * 
 * @author Daniel Seita
 */
class NormalHiddenMarkovModel {
    
    private int numStates;
    private int internalState;
    private double[][] transitions;                 // In log-space
    private NormalDistribution[] stateEmissions;    // NOT in log-space
    
    /** 
     * To initialize an HMM, call NornalHiddenMarkovModel(n, distr) where n is the number of states
     * and distr is a list of NormalDistribution objects that we can assign to the states. This HMM
     * class handles randomly initializing the transition probabilities and sampling observations.
     * 
     * @param numStates The number of states. We suggest 10 as a reasonable default.
     * @param stateEmissions An array of NormalDistributions, from the apache.commons.math3 library.
     */
    public NormalHiddenMarkovModel(int numStates, NormalDistribution[] stateEmissions) {
        this.numStates = numStates;
        this.stateEmissions = stateEmissions;
        this.internalState = 0;
        this.transitions = generateTransitionProbabilities(numStates);
        assert this.numStates == this.stateEmissions.length;
    }
    
    /** Returns the current state of the HMM, stored internally. */
    public int getCurrentState() {
        return this.internalState;
    }
    
    /** Resets the state back to zero, e.g., if we're in a loop and re-using this HMM */
    public void resetState() {
        this.internalState = 0;
    }
       
    /** Generates the next state according to the transition matrix (and exponentiates). */
    public int generateNextState() {
        double value = (new Random()).nextDouble();
        double sum = Math.exp(this.transitions[this.internalState][0]);
        for (int i = 0; i < numStates-1; i++) {
            if (value < sum) {
                this.internalState = i;
                return i;
            }
            sum += Math.exp(this.transitions[this.internalState][i+1]);
        }
        this.internalState = numStates-1;
        return numStates-1;
    }

    /** Samples an observation from our current state. */
    public double generateObservation() {
        return this.stateEmissions[internalState].sample();
    }
    
    /**
	 * Given the number of states, creates a transition matrix 'transitions' where transitions[k] is the
	 * vector of transition probabilities of the form a_{kx} where x is the index inside transitions[k],
	 * so k is the prior state and x is the current one. Right now it initializes values to be uniformly
	 * distributed between 0 and 1, and then normalizes by row. The transitions create an ergodic HMM.
	 * ALSO it then takes logs, so everything that uses it should be in log-space (or be aware of that).
	 * 
	 * @param numStates The number of states in the HMM, i.e., transitions.length and transitions[k].length.
	 * @return A random transition matrix, which could be used as the "A" matrix for an HMM.
	 */
	private double[][] generateTransitionProbabilities(int numStates) {
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
	            transitions[row][i] = Math.log(nonNormalizedValues[i] / sum);
	        }
	    }
	    return transitions;
	}
    
    /**
     * Given a list of observations generated (possibly NOT from this particular HMM), compute the forward
     * probability, i.e., the probability of generating that observation given this HMM. Note: we use log
     * probabilities, and also we only have the list of observations, so we need to sum over all possible
     * current states we could be in for the last observation.
     * 
     * Note: the emission probabilities are NOT assumed to be in log-space so take logs for those.
     * 
     * @param observations The list of observations, which are individually doubles.
     * @return The log of the forward probability.
     */
    public double getLogForwardProbability(List<Double> observations) {
        double result = Double.NEGATIVE_INFINITY;
        double[][] tempArray = new double[observations.size()][numStates]; // Holds cumulative (log) probabilities
        for (int state = 0; state < numStates; state++) {
            //double logProb = this.transitions[0][state] + this.emissions[state][observations.get(0)];
            double logProb = this.transitions[0][state] + Math.log(this.stateEmissions[state].density(observations.get(0)));
            assert !Double.isNaN(logProb) : "Problem: logProb = " + logProb;
            result = logAdd(result, logProb);
            tempArray[0][state] = logProb;
        }
        //assert Math.exp(result) >= 0.0 && Math.exp(result) <= 1.0 : "Problem: Math.exp(result) = " + Math.exp(result);
        if (observations.size() == 1) return result;
        for (int obs = 1; obs < observations.size(); obs++) {
            for (int j = 0; j < numStates; j++) {
                double logSum = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < numStates; i++) {
                    //double logComponent = tempArray[obs-1][i] + transitions[i][j] + emissions[j][observations.get(obs)];
                    double logComponent = tempArray[obs-1][i] + transitions[i][j] + 
                            Math.log(this.stateEmissions[j].density(observations.get(obs)));
                    assert !Double.isNaN(logComponent) : "logComponent = " + logComponent;
                    logSum = logAdd(logSum, logComponent);
                }
                tempArray[obs][j] = logSum;
            }
        }
        result = Double.NEGATIVE_INFINITY;
        for (int state = 0; state < numStates; state++) {
            result = logAdd(result, tempArray[observations.size()-1][state]);
        }
        return result;
    }

    /**
     * Given log(x) and log(y), returns log(x+y). This method is from Dan Klein, but the formula is
     * pretty commonplace.
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
