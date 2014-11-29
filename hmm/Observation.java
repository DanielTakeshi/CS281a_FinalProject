package hmm;

/**
 * 
 * 
 * @author Daniel Seita
 */
public class Observation {
    
    enum ObservationType {
        COMMUNICATION, TRUST, TRAVEL, MONEY
    }
    
    ObservationType type;
    
    public Observation(ObservationType type) {
        this.type = type;
    }
}
