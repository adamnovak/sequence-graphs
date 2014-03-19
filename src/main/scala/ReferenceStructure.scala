package edu.ucsc.genome

/**
 * Represents a Reference Structure: a phased or unphased sequence graph, with a
 * rule for mapping to it. Does not necessarily keep the SequenceGraph around.
 *
 * Internally, has two implementations: one for the bottom-level collection of
 * haplotypes, and one for a higher-level graph that merges nodes from the level
 * below.
 */
trait ReferenceStructure {
    /**
     * Map a string on both sides. Returns a sequence of Positions to which each
     * base maps, with None for bases that don't map. All the Position UUIDs
     * returned, which are (contig, base) pairs, will be specific to this level.
     */
    def map(pattern: String): Seq[Option[Position]]
    
    /**
     * Get the bottom-level FMDIndex, which upper-level ReferenceStructures need
     * in order to perform mapping operations efficiently.
     */
    def getIndex: FMDIndex
}

/**
 * A ReferenceStructure for phased sequence graphs (i.e. string haplotypes).
 * Backed by an FMDIndex.
 *
 * Takes a string basename for the index to load, as created with RLCSABuilder.
 * Positions are all on contigs in that index.
 */
class StringReferenceStructure(basename: String) extends ReferenceStructure {
    
    // Make our FMDIndex
    val index = new FMDIndex(basename)

    // Map with our index
    def map(pattern: String): Seq[Option[Position]] = index.map(pattern)

    // Expose our index to the level above us.
    def getIndex = index
}

/**
 * A ReferenceStructure which has some things in a lower-level
 * ReferenceStructure collapsed together. All positions in this
 * ReferenceStructure get unique IDs: Positions with base incrementing in the
 * order added and contig equal to the contig specified for the structure. TODO:
 * Re-design this class to allow multiple contigs inside it while still
 * providing shorthand to add new bases. Also re-design it to merge based on the
 * previous level rather than the lowest. And to efficiently persist to disk.
 */
class CollapsedReferenceStructure(base: ReferenceStructure, contig: String) 
    extends ReferenceStructure {

    // We have an IntervalMap of intervals we merge to what we merge them into.
    val intervals = new IntervalMap[Position]
    
    // We keep track of the next base ID that is free
    var nextId: Long = 1
    
    /**
     * Produce the next free base number.
     */
    def nextBase: Long = {
        // Keep the last value
        val toReturn = nextId
        // Increment
        nextId += 1
        // Return the value we kept
        toReturn
    }
    
    // We map using the ranges thing.
    def map(pattern: String): Seq[Option[Position]] = {
        // TODO: We're doing our own disambiguation here. Unify with the
        // FMDIndex disambiguate.
        
        println("Mapping with intervals:")
        println(intervals.mkString("\n"))
        
        // Map the bases in the string to optional range numbers, on the left
        val leftMappings = base.getIndex.leftMap(intervals.rangeVector, pattern)
        
        // Also map on the right, to completely different range numbers
        val rightMappings = base.getIndex.rightMap(intervals.rangeVector,
            pattern)
            
        // Convert to Positions by looking up in the interval value array.
        val leftPositions = leftMappings.map {
            case -1 => None
            // TODO: Change valueArray from a Seq to a map from long to
            // whatever, so it can hold a whole genome's worth of stuff.
            case rangeNumber => intervals.valueArray(rangeNumber.toInt)
        }
           
        val rightPositions = rightMappings.map {
            case -1 => None
            case rangeNumber => intervals.valueArray(rangeNumber.toInt)
        }
        
        println("Left positions:")
        println((leftMappings, leftPositions).zipped.mkString("\n"))
        println("Right positions:")
        println((rightMappings, rightPositions).zipped.mkString("\n"))
        println("Disambiguating...")
        
        // Disambiguate the alternative mappings, producing a sequence of final
        // mappings. TODO: Don't be calling down into a completely different
        // object.
        (leftPositions, rightPositions).zipped.map(getIndex.disambiguate(_, _))
    }
    
    // We return our base's index as our own.
    def getIndex = base.getIndex
    
    /**
     * Add a new Position at this level, created by merging the given sequence
     * of bottom-level positions. A given bottom-level position may be used only
     * once.
     */
    def addPosition(newPosition: Position, bottomPositions: Seq[Position]) = {
        bottomPositions.foreach { bottomPosition =>
            // Find the BWT index for this Position
            val bwtIndex = getIndex.inverseLocate(bottomPosition)
            
            // Store this new position in there, collecting with adjacent BWT
            // things if possible.
            intervals(bwtIndex, bwtIndex) = newPosition
        }
    }
    
    /**
     * Shorthand way to merge both Sides of a pair of positions into a new
     * position.
     */
    def merge(contig1: String, base1: Long, contig2: String, base2: Long) = {
        
        // Assign the next free base number to the merged position. It gets used
        // for both left and right sides, which make sense.
        val base = nextBase
        
        for(face <- Seq(Face.LEFT, Face.RIGHT)) {
            // Merge this face
            
            // What Position will we produce? It's in our structure's contig,
            // with the next free base number. TODO: We're inverting the face
            // because we need to store the positions correspondign to the left
            // faces of the things we merge (left-mapping semantics) as the
            // right faces of the things we merge them into (because of the
            // range-mapping functions' low-level right-mapping semantics.)
            val newPosition = new Position(contig, base, !face)
            
            // What positions is it based on?
            val toMerge = Seq(new Position(contig1, base1, face), 
                new Position(contig2, base2, face))
                
            // Actually add the position merging them together.
            addPosition(newPosition, toMerge)
        }
    }
    
    /**
     * Shorthand way of just taking a position from the bottom-level graph without
     * merging it with anything.
     */
    def pass(passContig: String, passBase: Long) = {
        // Assign the next free base number to the merged position. It gets used
        // for both left and right sides, which make sense.
        val base = nextBase
    
        for(face <- Seq(Face.LEFT, Face.RIGHT)) {
            // Pass this face
            
            // What Position will we re-name this to? TODO: We're inverting the
            // face because we need to store the positions correspondign to the
            // left faces of the things we pass (left-mapping semantics) as the
            // right faces of the things we pass them up as (because of the
            // range-mapping functions' low-level right-mapping semantics.)
            val position = new Position(contig, base, !face)
                
            // Actually add the position
            addPosition(position, Seq(new Position(passContig, passBase, face)))
        }
    }
}
