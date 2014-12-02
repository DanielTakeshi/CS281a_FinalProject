package hmm;

/**
 * 
 * 
 * @author Daniel Seita
 */
public class Entity {

    public enum EntityType {
        PERSON, PLACE, OBJECT
    }
    
    EntityType eType;
    
    public Entity(EntityType eType) {
        this.eType = eType;
    }

}
