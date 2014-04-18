package edu.ucsc.genome

import scala.collection.JavaConversions._
import scala.io.Source
import fi.helsinki.cs.rlcsa.{RangeVector, RangeEncoder, RangeVectorIterator}

// We use edu.ucsc.genome.Side objects to represent a base on a strand; face is
// the strand (LEFT for forward strand, RIGHT for reverse strand, under the idea
// that everything maps by upstream context). This is sort of backwards for
// suffix ranges, however.

/**
 * An index that allows efficient substring query and exact-match mapping to a
 * genome. The genome is indexed as a Heng Li-style FMD-Index: an FM-index
 * containing each sequence and its reverse complement. Internally, uses the
 * RLCSA FMD-Index extension via SWIG bindings.
 *
 * Can be serialized, but just keeps a reference to the filename of the on-disk
 * index, so it's not going to work very well across different machines unless
 * the index is in a consistent place.
 *
 * Supported operations:
 * - Basic FM-Index operations (count, locate)
 * - Mapping to genome Positions, or to specified BWT ranges, based on left or 
 *   right context.
 */
class FMDIndex(var basename: String) extends Serializable {
    import fi.helsinki.cs.rlcsa.{FMD, RLCSAUtil, MapAttemptResult, 
        MappingVector, Mapping, pair_type}
    
    if(basename == null) {
        throw new Exception("Cannot create FMDIndex with a null basename.")
    }
    
    // We keep a native FMD for the given basename, which we reload after every
    // serialization.
    @transient 
    lazy val fmd = new FMD(basename)
    
    // We load the contig file, keeping contig names by index, every time we
    // deserialize.
    @transient
    lazy val contigData: Array[(String, Long)] = Source.fromFile(basename +
        ".chrom.sizes").getLines.map { (line) =>
        
            // Split each line on the tab
            val parts = line.split("\t")
            // Make the second item a long
            (parts(0), parts(1).toLong)
        
        }.toArray
        
    // We also keep (index, length) by contig name. Also reloaded after
    // serialization. TODO: make this more efficient/use an on-disk database or
    // something.
    @transient
    lazy val contigInverse: Map[String, (Int, Long)] = 
        contigData.zipWithIndex.map {
            case ((contig: String, length: Long), index: Int) =>
                (contig, (index, length))
        }.toMap
     
    // Keep a bit vector so we can put in a Position ID and get out the number
    // of the contig it belongs to. There is a 1 at the start of the range
    // occupied by every contig. So the rank of an ID gives the contig it
    // belongs to. Strand doesn't matter here.
    @transient
    lazy val contigRangeVector: RangeVector = {
        // Make a new encoder with this arbitrary block size.
        val encoder = new RangeEncoder(32)
        
        // Keep track of the next ID to use. Assume the contigs take all the
        // first IDs.
        var nextID: Long = 0
        
        contigData.foreach { case (_, length) =>
            // For each contig
            
            // Put a 1 at the start of the contig
            encoder.addBit(nextID)
            
            // Push the next ID over by the length
            nextID += length
            
        }
        
        // Make the encoder be ready.
        encoder.flush()
        
        // Make the range vector
        new RangeVector(encoder, nextID + 1)
    }
    
    // Keep a persistent iterator for actual lookups. We keep the vector too
    // because it might manage to get GC'd if we don't, probably.
    @transient
    lazy val contigIdentifier: RangeVectorIterator = {
        new RangeVectorIterator(contigRangeVector)
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Metadata Operations
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Get all the contig names in the index, in order.
     */
    def contigs: Seq[String] = {
        contigData.map(_._1)
    }
        
    
    /**
     * Get the length of the contig with the given name.
     */
    def contigLength(contig: String): Long = {
        // We store (index, length) tuples by name.
        contigInverse(contig)._2
    }
    
    /**
     * Return the total length of all contigs.
     */
    def totalLength: Long = {
        // Sum up all the lengths
        contigData.map(_._2).sum
    }
    
    /**
     * Turn a 1-based contig name and base into a position ID.
     */
    def contigNameBaseToPosition(contig: String, base: Long): Long = {
        // What contig number is this?
        val contigNumber = contigInverse(contig)._1
        
        // What position ID beloings this far 1-based into it?
        contigNumberToPosition(contigNumber) + base - 1
    }
    
    /**
     * Get the contig number for a given Position number. Assumes that the given
     * position ID is within range of those assigned to contigs (i.e. less than
     * the total single-strand length of all contigs).
     */
    def positionToContigNumber(position: Long): Int = {
        // Just look up the rank of that position, and subtract for the 1 at the
        // start of the first contig.
        (contigIdentifier.rank(position) - 1).asInstanceOf[Int]
    }
    
    /**
     * Get the first ID used by any position on the contig with the given
     * number. This will be the the first base of the forward strand.
     */
    def contigNumberToPosition(contig: Int): Long = {
        // Ought to be the reverse of contigNumberForPosition.
        contigIdentifier.select(contig)
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // FM-Index Operations
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Return the number of occurrences of the given pattern in the index.
     */
    def count(pattern: String): Long = RLCSAUtil.length(fmd.count(pattern))
    
    /**
     * Return the Sides at which occurrences of the given pattern
     * start.
     */
    def locate(pattern: String): Iterator[Side] = {
        // Look up the range.
        val range = fmd.count(pattern)
        
        // Get a USIntArray of locations. We have to free it
        val locations = fmd.locate(range)
        
        // Make a list in our onw memory to copy to.
        var items: List[Side] = Nil
        
        for(i <- 0L until RLCSAUtil.length(range)) {
            // Get each result location. TODO: We can only go up to the max int
            // value in terms of number of locations.
            val location = RLCSAUtil.USIntArray_getitem(locations, i.toInt)
            
            // Get a text number and an offset
            val textAndOffset = fmd.getRelativePosition(location)
            
            // Make a Side and cons it onto the list.
            items = pairToSide(textAndOffset) :: items
        }
        
        // Free the C array of locations
        RLCSAUtil.delete_USIntArray(locations)
        
        // Give back an iterator for the list.
        // TODO: Do this iteration on demand.
        items.iterator
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Convert between Sides, text-offset pairs and BWT coordinates
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Convert a Side into a text number and offset in that text. If the
     * Side is a LEFT face, it will be on the forward strand; if it is a
     * RIGHT face, it will be on the reverse strand.
     */
    def sideToPair(side: Side): (Int, Long) = {
        // Look up the contig used in the side
        val contigNumber = positionToContigNumber(side.coordinate)
        // Look up its length
        val contigLength = contigData(contigNumber)._2
        // Get the base ID for the contig
        val contigStart = contigNumberToPosition(contigNumber)
        
        // What text index corresponds to the strand of the contig this Side
        // lives on? Forward texts are even, reverse ones are odd.
        val text = contigNumber * 2 + (if(side.face == Face.LEFT) 0 else 1)
        
        // What's the offset in the text (0-based)?
        val offset = if(side.face == Face.LEFT) {
            // We're on the forward strand, so just count up how far we have to
            // go out to get to that ID.
            side.coordinate - contigStart
        } else {
            // We're on the reverse strand, so we need to flip around and count
            // from the "start" of the reverse complement, and adjust to keep it
            // 0-based.
            contigLength - (contigStart - side.coordinate) - 1
        }
        
        // Return the text and offset
        (text, offset)
    }
    
    /**
     * Convert a Side to a BWT index in this index. If the Side is a
     * LEFT face, it will be on the forward strand; if it is a RIGHT face, it
     * will be on the reverse strand.
     */
    def sideToBWT(side: Side): Long = {
        // Turn the side into a text and offset.
        val (text, offset) = sideToPair(side)
        
        // Pack the two together in the right type to send down to C++
        val textAndOffset: pair_type = new pair_type(text, offset)
        
        // Turn the text and offset into an SA value in the combined all-text
        // coordinate space.
        val absolutePosition = fmd.getAbsolutePosition(textAndOffset)
        
        val relativePosition = fmd.getRelativePosition(absolutePosition)
        
        // Un-locate this SA value to find the index in the SA at which it
        // occurs.
        val SAIndex = fmd.inverseLocate(absolutePosition)
        
        // Convert to a BWT index and return
        SAIndex + fmd.getNumberOfSequences
    }
    
    /**
     * Convert a BWT index in this index to a Side. If the BWT index is on
     * the forward strand, it will be a Side with a LEFT face; if it is on
     * the reverse strand, it will be a Side with a RIGHT face.
     *
     * TODO: Consistent capitalization with function above?
     */
    def bwtToSide(bwt: Long): Side = {
        // Convert from BWT coordinates to SA coordinates, locate it, and then
        // report position as a (text, offset) pair.
        val textAndOffset = fmd.getRelativePosition(fmd.locate(bwt - 
            fmd.getNumberOfSequences));
        
        // Convert that pair to a Side
        pairToSide(textAndOffset)
    }
    
    /**
     * Convert a text and offset pair (which some C++-level mapping functions
     * use) to a Side. If the text is a forward- strand text, the Side's face
     * will be LEFT; otherwise it will be RIGHT.
     */
    def pairToSide(textAndOffset: pair_type): Side = {
        // Pull out the text and offset
        pairToSide(textAndOffset.getFirst.toInt, textAndOffset.getSecond)
    }
    
    /**
     * Convert a text number and offset length to a Side. If the text is a
     * forward- strand text, the Side's face will be LEFT; otherwise it will
     * be RIGHT.
     */
    def pairToSide(text: Int, offset: Long): Side = {
        // Work out how far we are into what strand of what contig.
        val contigNumber: Int = text / 2
        val contigStrand: Int = text % 2
        
        // Look up its length
        val contigLength = contigData(contigNumber)._2
        // Get the base ID for the contig
        val contigStart = contigNumberToPosition(contigNumber)
        
        // What face should we use for that strand?
        val face = contigStrand match {
            case 0 => Face.LEFT
            case 1 => Face.RIGHT
        }
        
        
        // What coordinate/Position ID should we use?
        val coordinate = contigStrand match {
            // On the forward strand, the first position is the first ID.
            case 0 => offset + contigStart
            // On the reverse strand we start at the end and go back. Correct to
            // stay 0-based.
            case 1 => contigLength - offset - 1
        }
        
        // Make a Side representing wherever we've decided this base really
        // goes.
        new Side(coordinate, face)
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Extracting sequence
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Return the base letter corresponding to the given Side.
     */
    def display(side: Side): Char = {
        // Find the text and offset for this side.
        val (text, offset) = sideToPair(side)
        
        // Make a single-item range.
        val rangeToGrab = new pair_type(offset, offset)
        
        println("About to display text %d offset %d".format(text, offset))
        
        // Grab the single-character array. This is really an array of shorts,
        // since Java bytes are signed and these are "unsigned chars", so a
        // short is needed to store the numerically-correct value.
        val array = fmd.display(text, rangeToGrab)
        
        
        // Get the (only) item in the array, and convert to to a Char
        RLCSAUtil.UCharArray_getitem(array, 0).asInstanceOf[Char]
        
        // Java should free the array I think...
        // TODO: Free the array?
    }
    
    
    ////////////////////////////////////////////////////////////////////////////
    // Left- and right-mapping single bases to Sides
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Return the side to which the given base in the given string uniquely
     * maps, or None if it does not uniquely map. Maps by context on the
     * specified face.
     *
     * Base is a 1-based index into the passed string.
     *
     * If the side's face is `Face.LEFT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.RIGHT`, it corresponds to
     * the reverse complement of the base it was mapped to.
     */
    def map(context: String, base: Int, face: Face): Option[Side] = {
        // Map a single base
        
        face match {
            case Face.RIGHT =>
                // Map the reverse complement on the other strand and flip
                // around.
                map(context.reverseComplement, context.size - (base - 1),
                    Face.LEFT) 
            case Face.LEFT =>
                // Map (0-based) and get a MapAttemptResult
                val mapping: MapAttemptResult = fmd.mapPosition(context,
                    base - 1)
                
                if(mapping.getIs_mapped) {
                    // Get the single actual mapped BWT position
                    val bwt = mapping.getPosition.getForward_start
                    
                    // Convert to a Side (for the very first character in
                    // the pattern)
                    val side = bwtToSide(bwt)
                    
                    // On the forward strand, we have to offset right by
                    // (mapping.getCharacters - 1), since we've gotten the Side
                    // of the leftmost character in the pattern and we really
                    // want that of this specific character. Unfortunately, on
                    // the reverse strand, we have to offset left by the same
                    // amount, since the pattern in that case is running
                    // backwards in genome coordinates.
                    val offset = side.face match {
                        case Face.LEFT => mapping.getCharacters - 1
                        case Face.RIGHT => -(mapping.getCharacters - 1)
                    }
                    
                    // Make a new Side, accounting for the offset from the
                    // left end of the pattern due to pattern length, and report
                    // that.
                    Some(new Side(side.coordinate + offset,
                        side.face))
                    
                } else {
                    None
                }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Left- and right-mapping entire strings to Sides per base.
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Map each Side in the given string by context on the given face,
     * returning a sequence of corresponding Sides, or None for bases that
     * don't map.
     *
     * If the Side's face is `Face.LEFT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.RIGHT`, it corresponds to
     * the reverse complement of the base it was mapped to.
     */
    def map(context: String, face: Face): Seq[Option[Side]] = {
        // Map a whole string
        
        face match {
            case Face.RIGHT =>
                // Map reverse complement on the other strand and flip around.
                map(context.reverseComplement, Face.LEFT).reverse
            case Face.LEFT =>
                // Do the actual mapping on this strand.
                val mappings = fmd.map(context)
                
                for {
                    i <- (0L until mappings.size)
                    mapping <- Some(mappings.get(i.toInt))
                } yield {
                    if(mapping.getIs_mapped) {
                        // Turn the (text, index) pair for this mapping into a
                        // Side. It's already been corrected for pattern
                        // length.
                        Some(pairToSide(mapping.getLocation))
                    } else {
                        // This didn't map, so the corresponding entry should be
                        // None
                        None
                    }
                }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Left- and right-mapping entire strings to range indices per base.
    ////////////////////////////////////////////////////////////////////////////

    /**
     * Map each side in the given string to a range starting with a 1 in the
     * given RangeVector, or to -1 if the FM-index search interval for the base
     * doesn't get contained in exactly 1 range. Uses the context on the given
     * face.
     *
     * The RangeVector specifies bi-ranges, where each range has a corresponding
     * reverse-complement range present. A range describes a right context.
     */
    def map(ranges: RangeVector, context: String, face: Face): Seq[Long] = {
        face match {
            case Face.LEFT =>
                // Map the reverse complement on the other strand and flip
                // around.
                map(ranges, context.reverseComplement, Face.RIGHT).reverse
            case Face.RIGHT =>
                // Do the actual mapping (on the right this time, since ranges
                // specify a downstream context).
                
                // FMD does all the work.
                val mappings = fmd.map(ranges, context)
                
                for {
                    i <- (0L until mappings.size)
                } yield {
                    // Convert vector to a Scala Seq
                    mappings.get(i.toInt)
                }
        }
    }
}

/**
 * Builds an RLCSA-format index at the given base name, merging things in if the
 * index already exists.
 *
 * Keeps the RLCSA files at the filenames that RLCSA derives from the given
 * basename, and the contig/node list at <basename>.chrom.sizes with <contig
 * name>\t<size> on each line, in the same order as the contigs appear in the
 * index. Note that the index contains each haplotype forwards and backwards,
 * but this file contains only one entry per haplotype, so e.g. index text 5 is
 * the reverse strand of the haplotype mentioned on line 3.
 *
 * Handles all the complexity of reverse-complementing FASTA records, turning
 * them into RLCSA null-terminated-string files, and maintaining a list of all
 * haplotypes by sequence number.
 */
class RLCSABuilder(basename: String) {
    import org.biojava3.core.sequence.io._
    import org.biojava3.core.sequence._
    import java.io._
    import java.nio.file._
    import scala.io._
    import scala.sys.process._
    import org.apache.commons.io._
    
    /**
     * Read all the sequences in the given FASTA file, index them and their
     * reverse complements, and merge that into the main index. The contig names
     * will be the FASTA IDs from the file.
     *
     * Note that this will place a couple copies of the FASTA data on /tmp.
     */
    def add(filename: String) = {
        // Get a temporary directory
        val scratchDirectory = Files.createTempDirectory("rlcsa")
        
        // Work out the name of the basename file that will store the haplotypes
        val haplotypes = scratchDirectory.resolve("haplotypes").toString
        // Open it up for writing
        val haplotypesWriter = new FileWriter(haplotypes)
        
        // Open the main index contig size list for appending
        val contigWriter = new FileWriter(basename + ".chrom.sizes", true)
    
        // Read all the FASTA records, lazily loading sequence
        val records: java.util.HashMap[String, DNASequence] = FastaReaderHelper
            .readFastaDNASequence(new File(filename), true)
        
        records.foreach { case (id, sequence) =>
            // Write each record to the haplotypes file
            // TODO: do we need to force capitalization here?
            haplotypesWriter.write(sequence.getSequenceAsString)
            // Terminated by a null
            haplotypesWriter.write("\0")
            
            // And then the reverse complement
            haplotypesWriter.write(sequence.getReverseComplement
                .getSequenceAsString)
            haplotypesWriter.write("\0")
            
            // Save the name and size in the sizes file
            contigWriter.write("%s\t%d\n".format(id, sequence.size))
        }
        
        haplotypesWriter.close
        contigWriter.close
        
        // TODO: Don't keep records in memory while doing the merge
        
        // Index the haplotypes file with build_rlcsa. Use some hardcoded number
        // of threads.
        Seq("build_rlcsa", haplotypes, "10").!
        
        // Now merge in the index
        merge(haplotypes)
        
        // Get rid of the temporary index files
        FileUtils.deleteDirectory(new File(scratchDirectory
            .toString))
        
    }
    
    /**
     * Merge in the index specified by otherBasename, or take it as our own
     * index if we don't have one yet. When this function returns, the
     * otherBasename files are no longer needed and can be deleted.
     *
     * DOES NOT handle merging chromosome size files; those are taken care of in
     * add().
     */
    private def merge(otherBasename: String) = {
        if(Files.exists(Paths.get(basename + ".rlcsa.array"))) {
            // We have an index already. Run a merge command.
            Seq("merge_rlcsa", basename, otherBasename, "10").!
            
        } else {
            // Take this index, renaming it to basename.whatever
            Files.copy(Paths.get(otherBasename + ".rlcsa.array"),
                Paths.get(basename + ".rlcsa.array"))
            Files.copy(Paths.get(otherBasename + ".rlcsa.parameters"),
                Paths.get(basename + ".rlcsa.parameters"))
            Files.copy(Paths.get(otherBasename + ".rlcsa.sa_samples"),
                Paths.get(basename + ".rlcsa.sa_samples"))
        }
    }
    
    /**
     * Print checksums of the built index to the console.
     */
    def checksum = {
        Seq("md5sum", basename + ".rlcsa.array").!
        Seq("md5sum", basename + ".rlcsa.parameters").!
        Seq("md5sum", basename + ".rlcsa.sa_samples").!
        Seq("md5sum", basename + ".chrom.sizes").!
    }
    
    /**
     * Produce an FMDIndex allowing you to query and map to the built index.
     * The FMDIndex will only work properly until the index is updated
     * again, at which point you should call this method again and get a new
     * one.
     */
    def getIndex: FMDIndex = new FMDIndex(basename)
    
}
