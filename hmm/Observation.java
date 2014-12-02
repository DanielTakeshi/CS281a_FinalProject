package hmm;

import java.util.*;

/**
 * An observation is the same thing as a transaction. Unlike normal HMMs, here an observation
 * consists of TWO states.
 * 
 * @author Daniel Seita
 */
public class Observation {
    
    public enum ObservationType {
        COMMUNICATION, TRUST, TRAVEL, MONEY
    }
    
    private ObservationType obsType;
    private State firstState;
    private State secondState;
    private List<Entity> people;
    private List<Entity> places;
    private List<Entity> objects;

    public Observation(State firstState, State secondState) {
        this.firstState = firstState;
        this.secondState = secondState;
        this.people = new ArrayList<Entity>();
        this.places = new ArrayList<Entity>();
        this.objects = new ArrayList<Entity>();
    }
}
