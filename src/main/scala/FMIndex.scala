package edu.ucsc.genome

import scala.collection.JavaConversions._
import scala.io.Source

// We use edu.ucsc.genome.Position objects to represent our positions. Contig is
// the actual haplotype name, base is the index from the start of the forward
// strand, and face is the strand (LEFT for forward strand, RIGHT for reverse
// strand, since everything maps by upstream context).

/**
 * A generic interface for implementations of an FMIndex and its associated
 * operations.
 * 
 */
trait FMIndex {

    /**
     * Return the number of occurrences of the given pattern in the index.
     */
    def count(pattern: String): Long

    /**
     * Return the positions at which occurrences of the given pattern
     * start.
     */
    def locate(pattern: String): Iterator[Position]
    
}

/**
 * An FMIndex that supports some more useful operations.
 */
trait FancyFMIndex extends FMIndex {
    
    /**
     * Return the minimum number of characters from the iterator required to get
     * exactly one unique match, or None if no unique match exists. Interprets
     * the iterator as specifying characters in reverse string order (so
     * characters will be pulled off and prepended to the start of the search
     * pattern).
     */
    def minUnique(characters: Iterator[Char]): Option[Long]
    
    /**
     * Return the position to which the given base in the given string uniquely
     * maps, or None if it does not uniquely map. Maps by upstream context (so
     * base 1, for example, probably won't map).
     *
     * Uses 1-based indexing for position and for base index.
     */
    def map(context: String, base: Int): Option[Position]

}

/**
 * An FMIndex that implements its fancier mapping features based on repeated
 * count operations rather than more efficient customized backwards search.
 *
 * Mapping on such an index will be O(n^2) instead of O(n) as it should be.
 */
trait BruteForceFMIndex extends FancyFMIndex {

    def minUnique(characters: Iterator[Char]): Option[Long] = {
        var pattern = ""
        
        while(count(pattern) > 1 && characters.hasNext) {
            // We haven't found a unique match, and we have more context. Use
            // it.
            pattern = characters.next + pattern
        }
        
        if(count(pattern) == 1) {
            // We found a unique match
            Some(pattern.size)    
        } else {
            // We found no match, or an ambiguous match and ran out of context.
            None
        }
        
    }
    
    def map(context: String, base: Int): Option[Position] = {
        // Map based on upstream context
        
        // Grab what's upstream. TODO: will this copy a whole chromosome?
        val upstreamContext = context.substring(0, base).reverse.iterator
        
        minUnique(upstreamContext) match {
            case Some(neededCharacters) => 
                // We can map with a certain amoutn of context. Grab that
                // context and map.
                
                // Grab the Position of the first base in the upstream context,
                // if we lay the last base down where it maps.
                val mappingPosition = locate(context.substring(
                    base - neededCharacters.toInt, base)).next
                
                // Fix it up to be the Position of the base we are actually
                // mapping (last one in the substring, since we map by upstream
                // context)
                mappingPosition.face match {
                    case Face.LEFT =>
                        // Go downstream
                        mappingPosition.base += neededCharacters
                    case Face.RIGHT =>
                        // Go upstream
                        mappingPosition.base -= neededCharacters
                }
                
                // Return the fixed up position to say we mapped.
                Some(mappingPosition)
            case None =>
                // We either ran out of context or found nothing matching.
                None
        }
    }
    
}

/**
 * An FMIndex that calls out to the "rlcsa_grep" CLI command to perform its
 * operations. Quite slow but requires only that rlcsa_grep be on the user's
 * PATH in order to work.
 *
 * Takes an RLCSA basename, which was the original name of the indexed file
 * (which may no longer really exist).
 */
class RLCSAGrepFMIndex(basename: String) extends BruteForceFMIndex {
    import scala.sys.process._

    // Load the contig size TSV, into an array of (name, size) by index.
    val chromSizes = Source.fromFile(basename + ".chrom.sizes").getLines.map { 
        (line) =>
        
        val parts = line.split("\t")
        // Parse out contig names and sizes.
        (parts(0), parts(1).toLong) 
    }.toArray

    def count(pattern: String): Long = {
        // Just run the command to count occurrences and intify its only output
        // line.
        Seq("rlcsa_grep", "-t", pattern, basename).!!.replace("\n", "").toLong
    }
    
    def locate(pattern: String): Iterator[Position] = {
        // Run the command to get all the occurrence positions (as (2 * contig
        // number + orientation), coordinate pairs) and asynchronously map its
        // resulting lines to longs. TODO: Use the -r option and look up contig
        // name and strand.
        Seq("rlcsa_grep", "-r", pattern, basename).lines.map { (line) =>
            // Split into text index and position
            val parts = line.replace("\n", "").split(", ")
            val textIndex = parts(0).toInt
            val position = parts(1).toLong
            
            // Compute the haplotype for the text index
            val contig: Int = textIndex / 2
            
            // Compute the orientation (0 for forward, 1 for reverse) that this
            // text means on that haplotype.
            val orientation = textIndex % 2
            
            // Look up the contig name and length
            val (contigName, contigLength) = chromSizes(contig)
            
            // Work out what the position on the forward strand is. This is the
            // first base in the pattern.
            val contigIndex = if(orientation == 0) {
                position
            } else {
                // Add 1 so that if we were at position contigLength we are now
                // at position 1.
                contigLength - position + 1
            }
            
            // Make a new Position talking about the right base of the right
            // strand of the right contig.
            new Position(contigName, contigIndex, 
                if(orientation == 0) Face.LEFT else Face.RIGHT)
        }.iterator
    }
}

/**
 * A class that uses the RLCSA FMD-Index extension via SWIG bindings.
 */
class FMDIndex(basename: String) extends BruteForceFMIndex {
    import fi.helsinki.cs.rlcsa.{FMD, RLCSAUtil, MapAttemptResult, 
        MappingVector, Mapping}
    
    // We keep a native FMD for the given basename
    val fmd = new FMD(basename)
    
    // We load the contig file
    val contigData: Array[(String, Long)] = Source.fromFile(basename +
        ".chrom.sizes").getLines.map { (line) =>
        
            // Split each line on the tab
            val parts = line.split("\t")
            // Make the second item a long
            (parts(0), parts(1).toLong)
        
        }.toArray
        
    // We keep the contig names
    
    /**
     * Get the Position for the last character in a pattern found at a given
     * offset in a given text, when the pattern has the given length. Returns
     * the Position for the start of that offset in that text (so left side of
     * base 1 for forward-strand offset 0, and right side of the last base for
     * reverse- strand offset 0).
     */
    def getPosition(text: Int, offset: Long, patternLength: Int): Position = {
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
        
        // What position should we use on that strand?
        val position = contigStrand match {
            // Forward strand we just go to the end of the pattern
            case 0 => offset + patternLength
            // On the reverse strand we start at the end and go back, but when
            // patternLength is 1 and offset is 0 we must return contigLength,
            // so we need the + 1.
            case 1 => contigLength - (offset + patternLength) + 1
        }
        
        // Make a Position representing wherever we've decided this base really
        // goes.
        new Position(contigName, position, face)
    }
    
    def count(pattern: String): Long = RLCSAUtil.length(fmd.count(pattern))
    
    def locate(pattern: String): Iterator[Position] = {
        // Look up the range.
        val range = fmd.count(pattern)
        
        // Get a USIntArray of locations. We have to free it
        val locations = fmd.locate(range)
        
        // Make a list in our onw memory to copy to.
        var items: List[Position] = Nil
        
        for(i <- 0L until RLCSAUtil.length(range)) {
            // Get each result location. TODO: We can only go up to the max int
            // value here.
            val location = RLCSAUtil.USIntArray_getitem(locations, i.toInt)
            
            // Get a text number and an offset
            val textAndOffset = fmd.getRelativePosition(location)
            
            // Cons it onto the list.
            items = getPosition(textAndOffset.getFirst.toInt,
                textAndOffset.getSecond, pattern.size) :: items
        }
        
        // Free the C array of locations
        RLCSAUtil.delete_USIntArray(locations)
        
        // Give back an iterator for the list.
        // TODO: Do this iteration on demand.
        return items.iterator
    }
    
    override def map(context: String, base: Int): Option[Position] = {
        // Map a single base
        
        // Map (0-based) and get a MapAttemptResult
        val mapping: MapAttemptResult = fmd.mapPosition(context, base - 1)
        
        if(mapping.getIs_mapped) {
            // Get the single actual mapped position
            val position = mapping.getPosition.getForward_start
            
            // Get a pointer to the same value
            val positionPointer = RLCSAUtil.copy_USIntPointer(position)
            
            // Convert to an SA coordinate
            fmd.convertToSAIndex(positionPointer);
            
            // Locate it, and then report position as a (text, offset) pair.
            // This will give us the position of the first base in the pattern,
            // which lets us infer the position of the last base in the pattern.
            val textAndOffset = fmd.getRelativePosition(fmd.locate(
                RLCSAUtil.USIntPointer_value(positionPointer)));
                
            // Free the pointer's memory
            RLCSAUtil.delete_USIntPointer(positionPointer);
                
            // Make a Position, accounting for the offset from the left end of
            // the pattern due to pattern length, and report that.
            Some(getPosition(textAndOffset.getFirst.toInt,
                textAndOffset.getSecond, mapping.getCharacters.toInt))
            
        } else {
            None
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
     * Produce an FMIndex allowing you to query the built index. The FMIndex
     * will only work properly until the index is updated again, at which point
     * you should call this method again and get a new one.
     */
    def getIndex: FMIndex = new RLCSAGrepFMIndex(basename)
    
}
