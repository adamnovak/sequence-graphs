package edu.ucsc.genome
import scala.collection.immutable.HashMap
import scala.collection.mutable.{ArrayBuilder, ArrayBuffer}
import org.ga4gh.{FMDUtil, BitVector, BitVectorIterator, FMDIndex, Mapping}
import scala.collection.JavaConversions._
import java.io.File
import java.nio.file._
import java.nio.{ByteBuffer, ByteOrder}
import java.nio.channels.FileChannel

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
     *
     * Ignores mappings on contexts shorter than the given minimum.
     */
    def map(context: String, minContext: Int = 0): Seq[Option[Side]] = {
        // Map on each side
        val leftMappings = mapFace(context, Face.LEFT, minContext)
        val rightMappings = mapFace(context, Face.RIGHT, minContext)
    
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
     *
     * Ignores mappings on contexts shorter than the given minimum.
     */
    def mapFace(context: String, face: Face, minContext: Int = 0):
        Seq[Option[Side]]
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
    def mapFace(pattern: String, face: Face, minContext: Int = 0):
        Seq[Option[Side]] = {
        
        face match {
            case Face.RIGHT =>
                // Do right-mapping as left-mapping flipped around.
                mapFace(pattern.reverseComplement, Face.LEFT, 
                    minContext).reverse
            case Face.LEFT =>
                // Get the MappingVector for mapping to all genomes with the
                // given minimum context.
                val mappings = getIndex.map(pattern, -1, minContext)
                
                // Make an ArrayBuffer of all the mappings, to which we can
                // efficiently append
                var mappingSeq: ArrayBuffer[Mapping] = new ArrayBuffer()
                
                for(i <- 0L until mappings.size()) {
                    mappingSeq += mappings.get(i.toInt)
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
 * Represents an array of Side objects that live in an on-disk file, as written
 * by createIndex (on an x86_64 system). The file format is a series of little-
 * endian 8 byte records, where the low bit represents the face (LEFT or RIGHT),
 * and the high 63 bits represent the ID. Unfortunately, this means in practice
 * we only get 63 bits of ID storage instead of 64.
 *
 */
class SideArray(filename: String) {
    // Open the file for random access.
    val file = FileChannel.open(Paths.get(filename), StandardOpenOption.READ)
    
    /**
     * Get the number of items in the array. This is just how many 8-byte
     * records fit in the file.
     */
    def length = file.size / 8
    
    /**
     * Dump the entire file in hex to the console, for debugging.
     */
    def dump = {
        // Make a 1-byte buffer to read into
        val buffer = ByteBuffer.allocate(1)
        
        (0L until length * 8).foreach { i =>
            // Grab each byte in turn
            file.read(buffer, i)
            
            // Rewind to start of buffer
            buffer.flip
            
            // Grab the byte as the first byte in the buffer
            val byteVal = buffer.get(0)
            
            // Dump the byte as hex
            println("%d: %02x".format(i, byteVal))
        }
    }
    
    /**
     * Load the record at the given index and return it as a Side.
     */
    def apply(index: Long): Side = {
        if(index > length) {
            // Don't let them look past the end.
            throw new Exception("Index %d beyond file length %d".format(index,
                length))
        }
        
        if(index < 0) {
            // Or before the beginning.
            throw new Exception("Index %d is negative".format(index))
        }
        
        // Make a buffer to read into
        val buffer = ByteBuffer.allocate(8)
        
        // Pop it into little endian mode
        buffer.order(ByteOrder.LITTLE_ENDIAN)
        
        // Read the bytes diurectly from the correct position in the file
        val bytesRead = file.read(buffer, index * 8)
        
        if(bytesRead != 8) {
            // Complain if we can't get all the bytes
            throw new Exception("Only managed to read " + bytesRead +
                " of 8 bytes")
        }
        
        // Rewind to start of buffer
        buffer.flip
        
        // Grab the record as the first long in the buffer.
        val record = buffer.getLong
        
        // Unpack the high 63 bits as the ID of the Side
        val coordinate = record >> 1
        // And the low bit as the face        
        val face = record & 1 match {
            case 0 => Face.LEFT
            case 1 => Face.RIGHT
        }
        
        // Make and return the Side
        new Side(coordinate, face)
        
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
        // Load the BitVector from it
        val toReturn = new BitVector(file)
        // Close the file. OS should free it.
        FMDUtil.fclose(file)
     
        // Send out the BitVector   
        toReturn
    }
    
    // Open (and wrap) the array of Sides that correspond to the ranges
    val sideArray = new SideArray(directory + "/mappings.bin")
        
    /**
     * Map the given string on the given side to all levels of the reference
     * structure. Ignore any mappings on less context than the specified minimum
     * context limit.
     */
    def mapFace(pattern: String, face: Face, minContext: Int = 0):
        Seq[Option[Side]] = {
        
        face match {
            case Face.LEFT =>
                // Try again on the right side.
                mapFace(pattern.reverseComplement, Face.RIGHT,
                    minContext).reverse
            case Face.RIGHT => 
                // Mapping to ranges is right-mapping.
                
                // Map to range numbers, or -1 for no mapping. This comes as a
                // SWIG- wrapped IntVector. Make sure to map to genome -1.
                val ranges = getIndex.map(rangeVector, pattern, -1, minContext)
                
                // Make an ArrayBuilder of all the mappings (which are Longs).
                // We use an ArrayBuilder instead of an ArrayBuffer since it's
                // better for primitive things like Longs. See
                // <http://stackoverflow.com/a/15839802/402891>
                val rangeSeq = new ArrayBuilder.ofLong
                
                for(i <- 0L until ranges.size()) {
                    rangeSeq += ranges.get(i.toInt)
                }
                
                // Convert to Sides in an Array and return.
                rangeSeq.result.map {
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
                            // We got a range that a Side is defined for.
                            // Retrieve and flip the Side.
                            Some(!(sideArray(range.toInt)))
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
        val iterator = new BitVectorIterator(rangeVector)
        // Look up the index and get the bit for every BWT position.
        for(i <- 0L until index.getBWTLength) yield {
            iterator.isSet(i)
        }
    }
}
