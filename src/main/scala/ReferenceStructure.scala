package edu.ucsc.genome
import scala.collection.immutable.HashMap

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
     * Given a Position on this level, get the BWT ranges corresponding to it.
     * At the bottom level, `Face.LEFT` positions are on the forward strand, and
     * `Face.RIGHT` positions are on the reverse. Ranges are inclusive.
     */
    def getRanges(position: Position): Seq[(Long, Long)]
    
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
 * Positions are all on contigs in that index.
 */
class StringReferenceStructure(index: FMDIndex) extends ReferenceStructure {
    
    /**
     * Make a new StringReferenceStructure from an index basename, describing an
     * FMD-index as would be created with RLCSABuilder.
     */
    def this(basename: String) = this(new FMDIndex(basename))

    // Map with our index
    def map(pattern: String, face: Face): Seq[Option[Position]] = {
        getIndex.map(pattern, face)
    }

    // Expose our index to the level above us.
    def getIndex = index
    
    def getRanges(position: Position): Seq[(Long, Long)] = {
        // Ranges for positions are just the BWT indices corresponding to them.
        val bwtPosition = index.positionToBWT(position)
        
        Seq((bwtPosition, bwtPosition))
    }
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
    
    // We also keep a map from parent Positions at this level to the lower-level
    // Positions they merge, which is necessary to get the ranges for a
    // Position. TODO: can we simplify this?
    var children: HashMap[Position, Seq[Position]] = HashMap.empty
    
    // We have an IDSource for base IDs
    val source = new IDSource(1)
    
    // We map using the range-based mapping mode on the index.
    def map(pattern: String, face: Face): Seq[Option[Position]] = {
        getIndex.map(intervals.rangeVector, pattern, face).map {
            // An range number of -1 means it didn't map 
            case -1 => None
            // Otherwise go get the Position for the range it mapped to (or None
            // if there's no position for that range). Make sure to flip it
            // around, to compensate for range mapping producing right-side
            // contexts on the forward strand instead of left-side ones.
            case range => intervals.valueArray(range.toInt).map(!_)
        }
    }
    
    // We return our base's index as our own.
    def getIndex = base.getIndex
    
    /**
     * Given a sequence of Positions from the previous level, create a new base
     * that merges all the given Positions into its left face, and the opposite
     * faces of all the given positions into its right face. Uses the specified
     * Position if given (must be a left Face), or creates a new Position
     * otherwise. Returns the Position of the left face of the base thus
     * created.
     */
    def addBase(toMerge: Seq[Position], target: Position = null): Position = {
        // First make a new position left face (i.e. forward orientation) if
        // needed.
        val newLeft = target match {
            case null => new Position(contig, source.id, Face.LEFT)
            case thing => thing
        }
        // And a flipped version for the right
        val newRight = !newLeft
        
        for(position <- toMerge; range <- base.getRanges(position)) {
            // Add in all the left side intervals pointing to our new left.
            intervals(range) = newLeft
        }
        
        for(position <- toMerge; range <- base.getRanges(!position)) {
            // Add in all the right side intervals pointing to our new right
            intervals(range) = newRight
        }
        
        // Add child lists
        children += ((newLeft, toMerge))
        children += ((newRight, toMerge.map(!_)))
        
        // Return our left-face Position
        newLeft
    }
    
    def getRanges(position: Position): Seq[(Long, Long)] = {
        // Just concatenate the ranges from all the clildren. TODO: Merging
        // ranges
        children(position).flatMap(base.getRanges(_))
    }
    
    /**
     * Add a new base using the given Position as its left side, and the
     * opposite face of that position's base as its right side.
     */
    def addBase(position: Position): Position = addBase(Seq(position))
    
    /**
     * Shorthand way to merge two bases from the lower level, in the same
     * orientation.
     */
    def merge(contig1: String, base1: Long, contig2: String, base2: Long) = {
        // Turn the two contig, base pairs into a sequence of Positions.
        addBase(Seq(new Position(contig1, base1, Face.LEFT),
            new Position(contig2, base2, Face.LEFT)))
    }
    
    /**
     * Shorthand way to pass a base from the lower level up to this level.
     */
    def pass(passContig: String, passBase: Long) = {
        // Make the contig, base pair into a position.
        addBase(new Position(passContig, passBase, Face.LEFT))
    }
}
