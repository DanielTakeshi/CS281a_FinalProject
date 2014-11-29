package hmm;

import java.util.*;

/**
 * I think I'm going to need to define the full sequence clearly. For now, have a set of Node objects
 * organized in a nice way ... I guess. Will need to think about this carefully
 * 
 * @author Daniel Seita
 */
public class ResearchDevelopmentHMM {
    
    public State root;

    /**
     * 
     */
    public ResearchDevelopmentHMM() {
        root = createRD();
    }

    /**
     * 
     * @return
     */
    public State createRD() {
        String[] stateNames = new String[] { "Initial military tactics and strategy discussions", 
                "Country X announces peaceful nuclear energy program", "Nuclear laboratory experiments",
                "Nuclear material acquisition for tests", "Enrichment pilot program detected",
                "Collect information on uranium enrichment", "Conduct implosion tests",
                "Acquisition of special test equipment", "Acquire produce material for high explosive tests",
                "Collect information of high explosives", "Purchase simulation software",
                "Request for declassified nuclear documents", "Country x solicits non-NPT countries for assistance",
                "Country X sends more students abroad in nuclear field", "Collect information on warhead design",
                "Scientists become active in nuclear community", "Leading scientists abroad recalled",
                "Start research and design activities"};
        State state = new State(stateNames[stateNames.length-1]);
        for (int i = stateNames.length-2; i >= 0; i--) {
            State newState = new State(stateNames[i]);
            newState.addChild(state, 0.5);
            state = newState;
        }
        return state;
    }
    
    /**
     * 
     */
    public void printRD() {
        root.print();
    }

}