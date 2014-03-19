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
     * Disambiguate a left mapping and a right mapping to produce an overall
     * mapping for a base, with left-mapping semantics (i.e. left face means
     * forward strand).
     */
    def disambiguate(leftMapping: Option[Position], 
        rightMapping: Option[Position]): Option[Position] = {
        
        (leftMapping, rightMapping) match {
            case (Some(leftPosition), Some(rightPosition)) =>
                // Only in this case do we have to check to ensure the mapping isn't
                // ambiguous.
                if(leftPosition.contig == rightPosition.contig && 
                    leftPosition.base == rightPosition.base && 
                    leftPosition.face == !rightPosition.face) {
                    
                    // They match (as opposite faces of the same base). Return the
                    // left one since we have left-mapping semantics.
                    Some(leftPosition)
                } else {
                    // They don't match. Ambiguous mappings don't count
                    None
                }
            case (Some(leftPosition), None) =>
                // Pass through the left mapping since we only have that.
                Some(leftPosition)
            case (None, Some(rightPosition)) => 
                // Right-only-mappings get faces flipped around to match our adopted
                // left-mapping semantics. Unary ! for Faces is defined in
                // package.scala.
                Some(new Position(rightPosition.contig, rightPosition.base, 
                !rightPosition.face))
            case (None, None) =>
                // It didn't map on either side.
                None
        }
    }
    
    /**
     * Map a string on both sides. Returns a sequence of Positions to which each
     * base maps, with None for bases that don't map. All the Position UUIDs
     * returned, which are (contig, base) pairs, will be specific to this level.
     *
     * If the position's face is `Face.LEFT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.RIGHT`, it corresponds to
     * the reverse complement of the base it was mapped to. This is a slight
     * semantic asymetry.
     * 
     * This is implemented by left-mapping and right-mapping the
     * base, and aggregating the results.
     */
    def map(context: String): Seq[Option[Position]] = {
        // Map on each side
        val leftMappings = map(context, Face.LEFT)
        val rightMappings = map(context, Face.RIGHT)
    
        // Zip them together and disambiguate each pair. Note that (a, b).zipped
        // is of a type that provides a map that takes binary functions, while
        // a.zip(b).map takes only unary functions.
        (leftMappings, rightMappings).zipped  map(disambiguate(_, _))
    }
    
    // Stuff to define in implementations
    
    /**
     * Get the bottom-level FMDIndex, which upper-level ReferenceStructures need
     * in order to perform mapping operations efficiently.
     */
    def getIndex: FMDIndex
    
    /**
     * Map all bases in the given string to Positions using the context on the
     * given face.
     */
    def map(context: String, face: Face): Seq[Option[Position]]
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
    def map(pattern: String, face: Face): Seq[Option[Position]] = {
        getIndex.map(pattern, face)
    }

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
    def map(pattern: String, face: Face): Seq[Option[Position]] = {
        getIndex.map(intervals.rangeVector, pattern, face).map {
            // An range number of -1 means it didn't map 
            case -1 => None
            // Otherwise go get the Position for the range it mapped to (or None
            // if there's no position for that range).
            case range => intervals.valueArray(range.toInt)
        }
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
            val bwtIndex = getIndex.positionToBWT(bottomPosition)
            
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
