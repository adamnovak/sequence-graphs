package edu.ucsc.genome
import scala.collection.immutable.HashMap
import org.ga4gh.{FMDUtil, RangeVector, RangeVectorIterator, FMDIndex, Mapping}
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
        face match {
            case Face.RIGHT =>
                // Do right-mapping as left-mapping flipped around.
                map(pattern.reverseComplement, Face.LEFT).reverse
            case Face.LEFT =>
                // Get the MappingVector
                val mappings = getIndex.map(pattern)
                
                // Make a Seq of all the mappings
                var mappingSeq: Seq[Mapping] = Seq()
                
                for(i <- 0L until mappings.size()) {
                    mappingSeq = mappingSeq :+ mappings.get(i.toInt)
                }
                
                // Grab the Sides that the mappings correspond to, or None.
                mappingSeq.map { mapping =>
                    if(mapping.getIs_mapped) {
                        // Go get the ID for this base, and match the
                        // appropriate Face. 
                        Some(new Side(index.getBaseID(mapping.getLocation), 
                            index.getStrand(mapping.getLocation) match { 
                                // Remember that the "true" strand is the
                                // reverse one.
                                case true => Face.RIGHT
                                case false => Face.LEFT
                            }))
                    } else {
                        // Didn't map anywhere.
                        None
                    }
                }
        }
    }

    // Expose our index to the level above us.
    def getIndex = index
    
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
        // We're about to use an API that doesn't really have error checking. So
        // first we make sure we can actually see this file.
        // TODO: check access rights and so forth.
        if(!(new File(directory + "/vector.bin").exists)) {
            throw new Exception("vector.bin file not found in %s"
                .format(directory))
        }
    
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
        face match {
            case Face.LEFT =>
                // Try again on the right side.
                map(pattern.reverseComplement, Face.RIGHT).reverse
            case Face.RIGHT => 
                // Mapping to ranges is right-mapping.
                
                // Map to range numbers, or -1 for no mapping. This comes as a
                // SWIG- wrapped IntVector.
                val ranges = getIndex.map(rangeVector, pattern)
                
                // Make a Seq of all the mappings
                var rangeSeq: Seq[Long] = Seq()
                
                for(i <- 0L until ranges.size()) {
                    rangeSeq = rangeSeq :+ ranges.get(i.toInt)
                }
                
                // Convert to Sides and return.
                rangeSeq.map {
                    // A range number of -1 means it didn't map 
                    case -1 => None
                    // Otherwise go get the Side for the range it mapped to (or
                    // None if there's no side for that range). Make sure to
                    // flip it around, to compensate for range mapping producing
                    // right-side contexts on the forward strand instead of
                    // left-side ones. Also remember that range indices are
                    // 1-based coming out of the FMD-index.
                    case range => 
                        if(range < sideArray.length + 1) {
                            // We got a range that a Side is defined for. Flip
                            // the Side.
                            Some(!(sideArray(range.toInt - 1)))
                        } else {
                            // Complain we're supposed to be mapping to a range
                            // that doesn't exist.
                            throw new Exception(
                                "Mapped to out-of-bounds range %d"
                                .format(range))
                        }
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
        for(i <- 0L until index.getBWTLength) yield {
            iterator.isSet(i)
        }
    }
}
