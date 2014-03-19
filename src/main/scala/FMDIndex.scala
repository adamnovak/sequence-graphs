package edu.ucsc.genome

import scala.collection.JavaConversions._
import scala.io.Source
import fi.helsinki.cs.rlcsa.RangeVector

// We use edu.ucsc.genome.Position objects to represent our positions. Contig is
// the actual haplotype name, base is the index from the start of the forward
// strand, and face is the strand (LEFT for forward strand, RIGHT for reverse
// strand, since everything maps by upstream context).

/**
 * An index that allows efficient substring query and exact-match mapping to a
 * genome. The genome is indexed as a Heng Li-style FMD-Index: an FM-index
 * containing each sequence and its reverse complement. Internally, uses the
 * RLCSA FMD-Index extension via SWIG bindings.
 *
 * Supported operations:
 * - Basic FM-Index operations (count, locate)
 * - Mapping to genome Positions, or to specified BWT ranges, based on left or 
 *   right context.
 */
class FMDIndex(basename: String) {
    import fi.helsinki.cs.rlcsa.{FMD, RLCSAUtil, MapAttemptResult, 
        MappingVector, Mapping, pair_type}
    
    if(basename == null) {
        throw new Exception("Cannot create FMDIndex with a null basename.")
    }
    
    // We keep a native FMD for the given basename
    val fmd = new FMD(basename)
    
    // We load the contig file, keeping contig names by index
    val contigData: Array[(String, Long)] = Source.fromFile(basename +
        ".chrom.sizes").getLines.map { (line) =>
        
            // Split each line on the tab
            val parts = line.split("\t")
            // Make the second item a long
            (parts(0), parts(1).toLong)
        
        }.toArray
        
    // We also keep (index, length) by contig name. TODO: make this more
    // efficient/use an on-disk database or something.
    val contigInverse: Map[String, (Int, Long)] = contigData.zipWithIndex.map {
        case ((contig: String, length: Long), index: Int) =>
            (contig, (index, length))
    }.toMap
    
    ////////////////////////////////////////////////////////////////////////////
    // FM-Index Operations
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Return the number of occurrences of the given pattern in the index.
     */
    def count(pattern: String): Long = RLCSAUtil.length(fmd.count(pattern))
    
    /**
     * Return the Positions at which occurrences of the given pattern
     * start.
     */
    def locate(pattern: String): Iterator[Position] = {
        // Look up the range.
        val range = fmd.count(pattern)
        
        // Get a USIntArray of locations. We have to free it
        val locations = fmd.locate(range)
        
        // Make a list in our onw memory to copy to.
        var items: List[Position] = Nil
        
        for(i <- 0L until RLCSAUtil.length(range)) {
            // Get each result location. TODO: We can only go up to the max int
            // value in terms of number of locations.
            val location = RLCSAUtil.USIntArray_getitem(locations, i.toInt)
            
            // Get a text number and an offset
            val textAndOffset = fmd.getRelativePosition(location)
            
            // Make a Position and cons it onto the list.
            items = pairToPosition(textAndOffset) :: items
        }
        
        // Free the C array of locations
        RLCSAUtil.delete_USIntArray(locations)
        
        // Give back an iterator for the list.
        // TODO: Do this iteration on demand.
        return items.iterator
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Convert between Positions, text-offset pairs and BWT coordinates
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Convert a Position to a BWT index in this index. If the Position is a
     * LEFT face, it will be on the forward strand; if it is a RIGHT face, it
     * will be on the reverse strand.
     */
    def positionToBWT(position: Position): Long = {
        // Look up the contig used in the Position
        val (contigNumber, contigLength) = contigInverse(position.contig)
        
        // What text index corresponds to the strand of the contig this Position
        // lives on? Forward texts are even, reverse ones are odd.
        val text = contigNumber * 2 + (if(position.face == Face.LEFT) 0 else 1)
        
        // What's the offset in the text?
        val offset = if(position.face == Face.LEFT) {
            // We're on the forward strand, so no change. Just convert from
            // 1-based Position coordinates to 0-based text coordinates.
            position.base - 1
        } else {
            // We're on the reverse strand, so we need to flip around and count
            // from the "start" of the reverse complement. This automatically
            // adjusts for the 1-based to 0-based conversion.
            contigLength - position.base
        }
        
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
     * Convert a BWT index in this index to a Position. If the BWT index is on
     * the forward strand, it will be a Position with a LEFT face; if it is on
     * the reverse strand, it will be a position with a RIGHT face.
     */
    def bwtToPosition(bwt: Long): Position = {
        // Convert from BWT coordinates to SA coordinates, locate it, and then
        // report position as a (text, offset) pair.
        val textAndOffset = fmd.getRelativePosition(fmd.locate(bwt - 
            fmd.getNumberOfSequences));
        
        // Convert that pair to a Position
        pairToPosition(textAndOffset)
    }
    
    /**
     * Convert a text and offset pair (which some C++-level mapping functions
     * use) to a Position. If the text is a forward- strand text, the position's
     * face will be LEFT; otherwise it will be RIGHT.
     */
    def pairToPosition(textAndOffset: pair_type): Position = {
        // Pull out the text and offset
        pairToPosition(textAndOffset.getFirst.toInt, textAndOffset.getSecond)
    }
    
    /**
     * Convert a text number and offset length to a Position. If the text is a
     * forward- strand text, the position's face will be LEFT; otherwise it will
     * be RIGHT.
     */
    def pairToPosition(text: Int, offset: Long): Position = {
        // Work out how far we are into what strand of what contig.
        val contigNumber: Int = text / 2
        val contigStrand: Int = text % 2
        
        // What contig name is that contig number?
        val contigName = contigData(contigNumber)._1
        // How long is that contig?
        val contigLength = contigData(contigNumber)._2
        
        // What face should we use for that strand?
        val face = contigStrand match {
            case 0 => Face.LEFT
            case 1 => Face.RIGHT
        }
        
        // What base number should we use on that strand?
        val base = contigStrand match {
            // On the forward strand, the first position is base 1
            case 0 => offset + 1
            // On the reverse strand we start at the end and go back.
            // Automatically corrects for 0/1-based-ness.
            case 1 => contigLength - offset
        }
        
        // Make a Position representing wherever we've decided this base really
        // goes.
        new Position(contigName, base, face)
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Left and right mapping to positions or BWT ranges
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Return the position to which the given base in the given string uniquely
     * maps, or None if it does not uniquely map. Maps by upstream context (so
     * base 1, for example, probably won't map).
     *
     * If the position's face is `Face.LEFT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.RIGHT`, it corresponds to
     * the reverse complement of the base it was mapped to.
     *
     * Uses 1-based indexing for position and for base index.
     */
    def leftMap(context: String, base: Int): Option[Position] = {
        // Map a single base
        
        // Map (0-based) and get a MapAttemptResult
        val mapping: MapAttemptResult = fmd.mapPosition(context, base - 1)
        
        if(mapping.getIs_mapped) {
            // Get the single actual mapped BWT position
            val bwt = mapping.getPosition.getForward_start
            
            // Convert to a Position (for the very first character in the
            // pattern)
            val position = bwtToPosition(bwt)
            
            // On the forward strand, we have to offset right by
            // (mapping.getCharacters - 1), since we've gotten the position of
            // the leftmost character in the pattern and we really want that of
            // this specific character. Unfortunately, on the reverse strand, we
            // have to offset left by the same amount, since the pattern in that
            // case is running backwards in genome coordinates.
            val offset = position.face match {
                case Face.LEFT => mapping.getCharacters - 1
                case Face.RIGHT => -(mapping.getCharacters - 1)
            }
            
            // Make a new Position, accounting for the offset from the left end
            // of the pattern due to pattern length, and report that.
            Some(new Position(position.contig, position.base + offset,
                position.face))
            
        } else {
            None
        }
        
    }
    
    /**
     * Map each position in the given string by upstream context, returning a
     * sequence of corresponding Positions, or None for bases that don't map.
     *
     * If the position's face is `Face.LEFT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.RIGHT`, it corresponds to
     * the reverse complement of the base it was mapped to.
     */
    def leftMap(context: String): Seq[Option[Position]] = {
        // Map a whole string
        
        // Do the actual mapping
        val mappings = fmd.map(context)
        
        for {
            i <- (0L until mappings.size)
            mapping <- Some(mappings.get(i.toInt))
        } yield {
            if(mapping.getIs_mapped) {
                // Turn the (text, index) pair for this mapping into a Position.
                // It's already been corrected for pattern length.
                Some(pairToPosition(mapping.getLocation))
            } else {
                // This didn't map, so the corresponding entry should be None
                None
            }
        }
    }
    
    /**
     * Map each position in the given string to a range starting with a 1 in the
     * given RangeVector, or -1 if the FM-index search interval for the base
     * doesn't get contained in exactly 1 range. Uses the left context.
     * 
     * This one comes pre-implemented in terms of right-mapping since the right-
     * sided one is more natural to implement on an FMD-index.
     */
    def leftMap(ranges: RangeVector, context: String): Seq[Long] = {
        // Flip the sequence around
        val reverseComplement = context.reverseComplement
        
        // Map it
        val mappings = rightMap(ranges, reverseComplement) 
        
        // Reverse the mappings (but leave the range numbers unchanged).
        mappings.reverse
    }
    
    /**
     * Return the position to which the given base in the given string uniquely
     * maps, or None if it does not uniquely map. Maps by downstream context.
     *
     * Uses 1-based indexing for position and for base index.
     *
     * If the position's face is `Face.RIGHT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.LEFT`, it corresponds to
     * the reverse complement of the base it was mapped to.
     *
     * By default, this is implemented by left-mapping the reverse complement.
     */
    def rightMap(context: String, base: Int): Option[Position] = {
        // Flip the sequence around
        val reverseComplement = context.reverseComplement
        
        // Left-map the corresponding base in the reverse complement.
        leftMap(reverseComplement, reverseComplement.size - (base - 1))
    }
    
    /**
     * Map each position in the given string by downstream context, returning a
     * sequence of corresponding Positions, or None for bases that don't map.
     *
     * If the position's face is `Face.RIGHT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.LEFT`, it corresponds to
     * the reverse complement of the base it was mapped to.
     *
     * By default, this is implemented by left-mapping the reverse complement.
     */
    def rightMap(context: String): Seq[Option[Position]] = {
        // Flip the sequence around
        val reverseComplement = context.reverseComplement
        
        // Map it
        val mappings = leftMap(reverseComplement) 
        
        // Reverse the mappings (but leave the Sides unchanged).
        mappings.reverse
    }
    
    /**
     * Map each position in the given string to a range starting with a 1 in the
     * given RangeVector, or -1 if the FM-index search interval for the base
     * doesn't get contained in exactly 1 range. Uses the right context.
     */
    def rightMap(ranges: RangeVector, context: String): Seq[Long] = {
        // Map a whole string to ranges with a range vector
        
        // FMD does all the work.
        val mappings = fmd.map(ranges, context)
        
        for {
            i <- (0L until mappings.size)
        } yield {
            // Convert vector to a Scala Seq
            mappings.get(i.toInt)
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Two-sided, ambiguity-proof mapping to positions only
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Disambiguate a left mapping and a right mapping to produce an overall
     * mapping for a base, with left-mapping semantics (i.e. left face means
     * forward strand).
     *
     * TODO: Make this something fancy with monads. TODO: Don't call down into
     * this from ReferenceStructure. Move this itno some utility place.
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
     * Map the given base in the given string on both sides. Returns the
     * position to which it maps, if it maps unambiguously, or None otherwise.
     *
     *
     * If the position's face is `Face.LEFT`, then the base mapped corresponds
     * to the base it was mapped to. If it is `Face.RIGHT`, it corresponds to
     * the reverse complement of the base it was mapped to. This is a slight
     * semantic asymetry.
     * 
     * By default, this is implemented by left-mapping and right-mapping the
     * base, and aggregating the results.
     */
    def map(context: String, base: Int): Option[Position] = {
        // Map on each side
        val leftMapping = leftMap(context, base)
        val rightMapping = rightMap(context, base)
        
        // Disambiguate the two mappings and return the result
        disambiguate(leftMapping, rightMapping)
    }
    
    /**
     * Map each position in the given string to a range starting with a 1 in the
     * given RangeVector, or -1 if the FM-index search interval for the base
     * doesn't get contained in exactly 1 range.
     *
     * By default, this is implemented by left-mapping and right-mapping the
     * string, and aggregating the results.
     */
    def map(context: String): Seq[Option[Position]] = {
        // Map on each side
        val leftMappings = leftMap(context)
        val rightMappings = rightMap(context)
    
        // Zip them together and disambiguate each pair. Note that (a, b).zipped
        // is of a type that provides a map that takes binary functions, while
        // a.zip(b).map takes only unary functions.
        (leftMappings, rightMappings).zipped  map(disambiguate(_, _))
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
     * Produce an FMDIndex allowing you to query and map to the built index.
     * The FMDIndex will only work properly until the index is updated
     * again, at which point you should call this method again and get a new
     * one.
     */
    def getIndex: FMDIndex = new FMDIndex(basename)
    
}
