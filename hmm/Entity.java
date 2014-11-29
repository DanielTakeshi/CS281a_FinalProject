package hmm;

/**
 * 
 * 
 * @author Daniel Seita
 */
public class Entity {

    enum EntityType {
        PERSON, PLACE, OBJECT
    }
    
    EntityType entityType;
    
    public Entity(EntityType entityType) {
        this.entityType = entityType;
    }

}
