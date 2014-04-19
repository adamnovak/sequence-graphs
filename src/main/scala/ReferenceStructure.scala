package edu.ucsc.genome
import scala.collection.immutable.HashMap
import fi.helsinki.cs.rlcsa.{FMDUtil, RangeVector, RangeVectorIterator}
import org.apache.avro.specific.{SpecificDatumReader, SpecificRecord}
import org.apache.avro.file.DataFileReader
import scala.collection.JavaConversions._
import java.io.File

/**
 * Represents a Reference Structure: a phased or unphased sequence graph, with a
 * rule for mapping to it. Does not necessarily keep the SequenceGraph around.
 *
 * Internally, has two implementations: one for the bottom-level collection of
 * haplotypes, and one for a higher-level graph that merges nodes from the level
 * below.
 *
 * Must be serializable, but re-loading the serialized version can depend on
 * references to external index files (like the FMD-index).
 */
trait ReferenceStructure extends Serializable {
    /**
     * Disambiguate a left mapping and a right mapping to produce an overall
     * mapping for a base, with left-mapping semantics (i.e. left face means
     * forward strand).
     */
    def disambiguate(leftMapping: Option[Side], 
        rightMapping: Option[Side]): Option[Side] = {
        
        (leftMapping, rightMapping) match {
            case (Some(leftSide), Some(rightSide)) =>
                // Only in this case do we have to check to ensure the mapping
                // isn't ambiguous.
                if(leftSide  == !rightSide) {
                    // They match (as opposite faces of the same base). Return
                    // the left one since we have left-mapping semantics.
                    Some(leftSide)
                } else {
                    // They don't match. Ambiguous mappings don't count
                    None
                }
            case (Some(leftSide), None) =>
                // Pass through the left mapping since we only have that.
                Some(leftSide)
            case (None, Some(rightSide)) => 
                // Right-only-mappings get faces flipped around to match our
                // adopted left-mapping semantics.
                Some(!rightSide)
            case (None, None) =>
                // It didn't map on either side.
                None
        }
    }
    
    /**
     * Map a string on both sides. Returns a sequence of Sides to which each
     * base maps, with None for bases that don't map. All the Side UUIDs
     * returned, which are position longs, will be specific to this level.
     *
     * If the side's face is `Face.LEFT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.RIGHT`, it corresponds to
     * the reverse complement of the base it was mapped to. This is a slight
     * semantic asymetry.
     * 
     * This is implemented by left-mapping and right-mapping the
     * base, and aggregating the results.
     */
    def map(context: String): Seq[Option[Side]] = {
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
     * Given a Side on this level, get the BWT ranges corresponding to it.
     * At the bottom level, `Face.LEFT` sides are on the forward strand, and
     * `Face.RIGHT` sides are on the reverse. Ranges are inclusive.
     */
    def getRanges(side: Side): Seq[(Long, Long)]
    
    /**
     * Map all bases in the given string to Sides using the context on the
     * given face.
     */
    def map(context: String, face: Face): Seq[Option[Side]]
}

/**
 * A ReferenceStructure for phased sequence graphs (i.e. string haplotypes).
 * Backed by an FMDIndex.
 *
 * Sides are all on bottom-level contigs in that index.
 */
class StringReferenceStructure(index: FMDIndex) extends ReferenceStructure {
    
    /**
     * Make a new StringReferenceStructure from an index basename, describing an
     * FMD-index as would be created with RLCSABuilder.
     */
    def this(basename: String) = this(new FMDIndex(basename))

    // Map with our index
    def map(pattern: String, face: Face): Seq[Option[Side]] = {
        getIndex.map(pattern, face)
    }

    // Expose our index to the level above us.
    def getIndex = index
    
    def getRanges(side: Side): Seq[(Long, Long)] = {
        // Ranges for sides are just the BWT indices corresponding to them.
        val bwtSide = index.sideToBWT(side)
        
        Seq((bwtSide, bwtSide))
    }
}

/**
 * A ReferenceStructure that has been built by the createIndex program and
 * loaded form disk. Internally keeps track of its bit vector of ranges and
 * array of Sides to which ranges map things.
 */
class MergedReferenceStructure(index: FMDIndex, directory: String)
    extends ReferenceStructure {
    
    // Load the range vector
    val rangeVector = {
        // First open a C FILE* with the minimal API that FMDUtil in RLCSA comes
        // with.
        val file = FMDUtil.fopen(directory + "/vector.bin", "r")
        // Load the RangeVector from it
        val toReturn = new RangeVector(file)
        // Close the file. OS should free it.
        FMDUtil.fclose(file)
     
        // Send out the RangeVector   
        toReturn
    }
    
    // Load the array of Sides that correspond to the ranges
    val sideArray: Array[Side] = {
        // We're loading a plain Avro file (no Parquet). See <http://avro.apache
        // .org/docs/current/gettingstartedjava.html#Deserializing>
    
        // Make a DatumReader because Avro thinks we might want to not have one
        // sometimes.
        val datumReader = new SpecificDatumReader[Side](classOf[Side])
        
        // Make a file reader to read the file with the datum reader.
        val fileReader = new DataFileReader[Side](
            new File(directory + "/mappings.avro"), datumReader)
        
        // Get an iterator over the file, and turn it into an array.
        fileReader.iterator.toArray
    }
        
    // We map using the range-based mapping mode on the index.
    def map(pattern: String, face: Face): Seq[Option[Side]] = {
        // Map to defined ranges in the range vector
        getIndex.map(rangeVector, pattern, face).map {
            // An range number of -1 means it didn't map 
            case -1 => None
            // Otherwise go get the Side for the range it mapped to (or None
            // if there's no side for that range). Make sure to flip it
            // around, to compensate for range mapping producing right-side
            // contexts on the forward strand instead of left-side ones.
            case range => 
                if(range < sideArray.length) {
                    // We got a range that a Side is defined for. Flip the Side.
                    Some(!(sideArray(range.toInt)))
                } else {
                    // We probably got the 1-past-the-end range? TODO: Did we
                    // break something?
                    None
                }
        }
    }
    
    def getIndex: FMDIndex = index
    
    /**
     * Dump the whole BitVector.
     */
    def bits: Seq[Boolean] = {
        val iterator = new RangeVectorIterator(rangeVector)
        // Look up the index and get the bit for every BWT position.
        for(i <- 0L until index.bwtRange._2) yield {
            iterator.isSet(i)
        }
    }
    
    /**
     * Get all the ranges in order, with their Sides.
     */
    def getRanges: Seq[((Long, Long), Side)] = {
        // Get an iterator
        val iterator = new RangeVectorIterator(rangeVector)
    
        for((side, i) <- sideArray.zipWithIndex)  yield {
            ((iterator.select(i), iterator.select(i + 1) - 1), side)
        }
    }
    
    def getRanges(side: Side): Seq[(Long, Long)] = {
        
        // Get an iterator
        val iterator = new RangeVectorIterator(rangeVector)
        
        for(
            (candidate, index) <- sideArray.zipWithIndex;
            if side == candidate
        ) yield {
            // This index is one of the ranges we want. Get its bounds.
            (iterator.select(index), iterator.select(index + 1) - 1)
        }
    }
}

/**
 * A ReferenceStructure which has some things in a lower-level
 * ReferenceStructure collapsed together. All sides in this ReferenceStructure
 * get unique IDs: Sides with base incrementing in the order added. TODO: Re-
 * design to merge based on the previous level rather than the lowest.
 */
class CollapsedReferenceStructure(base: ReferenceStructure) 
    extends ReferenceStructure {
    
    // We have an IntervalMap of intervals we merge to what we merge them into.
    val intervals = new IntervalMap[Side]
    
    // We also keep a map from parent Sides at this level to the lower-level
    // Sides they merge, which is necessary to get the ranges for a
    // Side. TODO: can we simplify this?
    var children: HashMap[Side, Seq[Side]] = HashMap.empty
    
    // We have an IDSource for base IDs. TODO: not safe across multiple levels.
    val source = new IDSource(getIndex.totalLength)
    
    // We map using the range-based mapping mode on the index.
    def map(pattern: String, face: Face): Seq[Option[Side]] = {
        getIndex.map(intervals.rangeVector, pattern, face).map {
            // An range number of -1 means it didn't map 
            case -1 => None
            // Otherwise go get the Side for the range it mapped to (or None
            // if there's no side for that range). Make sure to flip it
            // around, to compensate for range mapping producing right-side
            // contexts on the forward strand instead of left-side ones.
            case range => 
                if(range < intervals.valueArray.length) {
                    // We got a range that an interval is defined for.
                    intervals.valueArray(range.toInt).map(!_)
                } else {
                    // We probably got the 1-past-the-end range. Maybe our
                    // IntervalMap is empty?
                    None
                }
        }
    }
    
    // We return our base's index as our own.
    def getIndex = base.getIndex
    
    /**
     * Given a sequence of Sides from the previous level, create a new base
     * that merges all the given Sides into its left face, and the opposite
     * faces of all the given sides into its right face. Uses the specified
     * Side if given (must be a left Face), or creates a new Side
     * otherwise. Returns the Side of the left face of the base thus
     * created.
     */
    def addBase(toMerge: Seq[Side], target: Side = null): Side = {
        // First make a new side left face (i.e. forward orientation) if
        // needed.
        val newLeft = target match {
            case null =>
                // We need to allocate a new ID.
                new Side(source.id, Face.LEFT)
            case thing => thing
        }
        // And a flipped version for the right
        val newRight = !newLeft
        
        for(side <- toMerge; range <- base.getRanges(side)) {
            // Add in all the left side intervals pointing to our new left.
            intervals(range) = newLeft
        }
        
        for(side <- toMerge; range <- base.getRanges(!side)) {
            // Add in all the right side intervals pointing to our new right
            intervals(range) = newRight
        }
        
        // Add child lists
        children += ((newLeft, toMerge))
        children += ((newRight, toMerge.map(!_)))
        
        // Return our left-face Side
        newLeft
    }
    
    def getRanges(side: Side): Seq[(Long, Long)] = {
        // Just concatenate the ranges from all the clildren. TODO: Merging
        // ranges
        children(side).flatMap(base.getRanges(_))
    }
    
    /**
     * Add a new base using the given Side as its left side, and the
     * opposite face of that side's base as its right side.
     */
    def addBase(side: Side): Side = addBase(Seq(side))
    
    /**
     * Shorthand way to merge two bases from the lower level, in the same
     * orientation.
     */
    def merge(contig1: String, base1: Long, contig2: String, base2: Long) = {
        // Turn the two contig, base pairs into a sequence of Sides.
        
        // First turn them into position IDs.
        val pos1 = getIndex.contigNameBaseToPosition(contig1, base1)
        val pos2 = getIndex.contigNameBaseToPosition(contig2, base2)

        // Then make Sides and add the bases.        
        addBase(Seq(new Side(pos1, Face.LEFT), new Side(pos2, Face.LEFT)))
    }
    
    /**
     * Shorthand way to pass a base from the lower level up to this level.
     */
    def pass(passContig: String, passBase: Long) = {
        // Make the contig, base pair into a side.
        addBase(new Side(
            getIndex.contigNameBaseToPosition(passContig, passBase), Face.LEFT))
    }
}
