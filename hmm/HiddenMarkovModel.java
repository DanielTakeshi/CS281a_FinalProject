package hmm;

import java.util.*;
import utilities.*;

/**
 * Class for HMMs I'll use. Keep in mind that calculations are to be done in log space, so be sure to convert
 * to normal-space if it makes things easier.
 * 
 * @author Daniel Seita
 */
class HiddenMarkovModel {
    
    private int numStates;
    private int internalState;
    public double[][] transitions;
    public double[][] emissions;
    
    public HiddenMarkovModel(double[][] transitions, double[][] emissions) {
        this.numStates = transitions.length;
        this.internalState = 0;
        this.transitions = transitions;
        this.emissions = emissions;
    }
    
    /**
     * @return The current state of the HMM, stored internally
     */
    public int getCurrentState() {
        return this.internalState;
    }
       
    /**
     * Generates the next state according to the transition matrix. I'll exponentiate to make things easier.
     * 
     * @return The next state, probabilistically chosen
     */
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
    
    /**
     * Probabilistically generates a binary-valued observation depending on the state we're in. This is
     * easiest to implement when comparing with the exponentiated version of this.emissions[state][0].
     * 
     * @return The observation generated from the state, probabilistically chosen
     */
    public int generateObservation() {
        if ((new Random()).nextDouble() < Math.exp(this.emissions[this.internalState][0])) {
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
            result = DanielUtilities.logAdd(result, logProb);
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
                    logSum = DanielUtilities.logAdd(logSum, logComponent);
                }
                tempArray[obs][j] = logSum;
            }
        }
        result = Double.NEGATIVE_INFINITY;
        for (int state = 0; state < numStates; state++) {
            result = DanielUtilities.logAdd(result, tempArray[observations.size()-1][state]);
        }
        return result;
    }
}
