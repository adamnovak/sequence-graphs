package edu.ucsc.genome

/**
 * Represents a thing that can create long IDs sequentially. Mutable.
 */
class IDSource(idStart: Long = 0, idCount: Long = -1) extends Serializable {

    // Compute the actual ID limit, which may or many not be present, but if
    // present should be what the caller passed.
    val idLimit = idCount match {
        case -1 => None
        case other => Some(idStart + other)
    }

    // Keep track of the next ID available to allocate
    var nextID = idStart
    
    /**
     * Allocate and return an ID. If we have run out of IDs, raise an exception.
     */
    def id: Long = {
        // Grab the next ID value
        val toReturn = nextID
        // Say the next next ID is 1 after that one
        nextID += 1
        
        idLimit.map { limit =>
            if(toReturn >= limit) {
                // Complain we ran out of IDs
                throw new Exception(
                    "Source can't allocate %d from ID block of %d at %d".format(
                    toReturn, idCount, idStart))
            }
        }
        
        // Return the ID we generated
        toReturn
    }
    
    /**
     * Allocate a whole block of the given number of IDs. If we would
     * run out of IDs, raise an exception.
     *
     * Return the first ID allocated; the block is contiguous, extending up from
     * there.
     */
    def ids(count: Long): Long = {
        // What's the first ID we're going to use?
        val firstID = nextID
        
        // Skip ahead to after this whole block.
        nextID += count
        
        idLimit.map { limit =>
            // We have a limit
            if(nextID > limit) {
                // We've run out.
                throw new Exception(
                    "Source can't allocate %d to %d from ID block of %d at %d"
                    .format(firstID, nextID, idCount, idStart))
            }
        }
        
        // Return the first ID
        firstID
    }
    
    /**
     * Split off a new IDSource that can make the given number of IDs.
     */
    def split(ids: Long): IDSource = {
        idLimit.map { (limit) =>
            // TODO: check bounds here
            if(nextID + ids >= limit) {
                // Complain we might run out of IDs
                throw new Exception(
                    "Source can't reserve %d from ID block of %d at %d".format(
                    ids, idCount, idStart))
            }
        }
        
        // If we get here we won't hit a limit. Make a new source.
        val toReturn = new IDSource(nextID, ids)
        
        // Skip the IDs it gets
        nextID = nextID + ids
        
        // Return it
        toReturn
    }
    
    /**
     * Split off a certain number of ID sources that can each make a certain
     * number of IDs.
     */
    def split(ids: Long, count: Long): Seq[IDSource] = {
        for(i <- 0L until count) yield split(ids)
    }
}
