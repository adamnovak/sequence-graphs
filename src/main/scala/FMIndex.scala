package edu.ucsc.genome

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
     * base 0, for example, probably won't map).
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
        val upstreamContext = context.substring(0, base + 1).reverse.iterator
        
        minUnique(upstreamContext) match {
            case Some(neededCharacters) => 
                // We can map with a certain amoutn of context. Grab that
                // context and map.
                
                // TODO: Roll all these operations together at the C level so we
                // can just backwards search until we have a unique match or
                // fail.
                // TODO: Make this work with longs.
                Some(locate(context.substring(base - neededCharacters.toInt + 1, 
                    base + 1)).next)
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

    def count(pattern: String): Long = {
        // Just run the command to count occurrences and intify its only output
        // line.
        Seq("rlcsa_grep", "-t", pattern, basename).!!.replace("\n", "").toLong
    }
    
    def locate(pattern: String): Iterator[Position] = {
        // Run the command to get all the occurrence positions (in global
        // concatenated coordinates), and asynchronously map its resulting lines
        // to longs.
        Seq("rlcsa_grep", "-s", pattern, basename).lines.map { (line) =>
            new Position("index", line.replace("\n", "").toLong, Face.LEFT)
        }.iterator
    }
}
