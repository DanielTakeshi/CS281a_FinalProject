package hmm;

import java.util.*;

/**
 * 
 * 
 * @author Daniel Seita
 */
class State {
    
    private String title;
    private List<State> children;
    private List<Double> childProbabilities;
    private double selfLoopProbability;
    private boolean isTerminal;
    
    public State(String title) {
        this.title = title;
        this.children = new ArrayList<State>();
        this.childProbabilities = new ArrayList<Double>();
        this.selfLoopProbability = 1.0;
        this.isTerminal = true;
    }

    public int numChildren() {
        return this.children.size();
    }
    
    public double getSelfLoopProb() {
        return this.selfLoopProbability;
    }
    
    /**
     * @return The list of children, not including itself.
     */
    public List<State> getChildren() {
        return this.children;
    }
    
    /**
     * @return The list of probabilities, not including the probability of a loop-back to itself
     */
    public List<Double> getProbabilities() {
        return this.childProbabilities;
    }

    public boolean isEndState() {
        return this.isTerminal;
    }
    
    /**
     * Adds a child to the list of children. Here's how the normalization works:
     * we normalize all other probabilities NOT including this current one. Then
     * we set the new child to have transition probability of "transitionProbability"
     * That way, when we add a state, the probability we put there IS the
     * probability we'll get. Notice that this will also make "isTerminal"
     * false, if it is not done so already.
     * 
     * @param s
     * @param transitionProbability
     */
    public void addChild(State s, double transitionProbability) {
        if (transitionProbability > 1.0 || transitionProbability < 0.0) {
            System.out.println("ERROR! Assigning a transition probability of " + transitionProbability + ".");
            System.exit(0);
        }
        double scalingFactor = 1 - transitionProbability;
        this.selfLoopProbability = this.selfLoopProbability * scalingFactor;
        for (int i = 0; i < this.childProbabilities.size(); i++) {
            this.childProbabilities.set(i, this.childProbabilities.get(i) * scalingFactor);
        }
        this.children.add(s);
        this.childProbabilities.add(transitionProbability);
        this.isTerminal = false;
    }
    
    @Override
    public String toString() {
        String result = "{{{ Tile = " + this.title + ", Self-loop prob = " + this.selfLoopProbability
                + ", Children probs = " + this.getProbabilities() + " }}}";
        return result;
    }

    /**
     * Use this to pretty-print the tree on the command line. It recursively calls the other print() method.
     */
    public void print() {
        print("", true);
    }

    /**
     * 
     * @param prefix
     * @param isTail
     */
    private void print(String prefix, boolean isTail) {
        System.out.println(prefix + (isTail ? "\\----------> " : "|----------> ") + this.title);
        for (int i = 0; i < children.size() - 1; i++) {
            children.get(i).print(prefix + (isTail ? "    " : "|    "), false);
        }
        if (children.size() > 0) {
            children.get(children.size() - 1).print(prefix + (isTail ?"    " : "|    "), true);
        }
    }
}