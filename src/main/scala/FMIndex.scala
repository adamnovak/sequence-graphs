package edu.ucsc.genome

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
    def locate(pattern: String): Iterator[Long]
    
}

/**
 * An FMIndex that supports some more useful operations.
 */
trait FancyFMIndex extends FMIndex {
    
    /**
     * Return the minimum number of characters from the iterator required to get
     * exactly one unique match, or None if no unique match exists. Interprets
     * the iterator as specifying characters in string order (so characters will
     * be pulled off and appended to the end of the search pattern).
     */
    def minUnique(characters: Iterator[Char]): Option[Long]

}

/**
 * An FMIndex that calls out to the "rlcsa_grep" CLI command to perform its
 * operations. Quite slow but requires only that rlcsa_grep be on the user's
 * PATH in order to work.
 *
 * Takes an RLCSA basename, which was the original name of the indexed file
 * (which may no longer really exist).
 */
class RLCSAGrepFMIndex(basename: String) extends FancyFMIndex {
    import scala.sys.process._

    def count(pattern: String): Long = {
        // Just run the command to count occurrences and intify its only output
        // line.
        Seq("rlcsa_grep", "-t", pattern, basename).!!.replace("\n", "").toLong
    }
    
    def locate(pattern: String): Iterator[Long] = {
        // Run the command to get all the occurrence positions (in global
        // concatenated coordinates), and asynchronously map its resulting lines
        // to longs.
        Seq("rlcsa_grep", "-s", pattern, basename).lines.map {
            _.replace("\n", "").toLong
        }.iterator
    }
    
    def minUnique(characters: Iterator[Char]): Option[Long] = {
        var pattern = ""
        
        while(count(pattern) > 1 && characters.hasNext) {
            // We haven't found a unique match, and we have more context. Use
            // it.
            pattern += characters.next
        }
        
        if(count(pattern) == 1) {
            // We found a unique match
            Some(pattern.size)    
        } else {
            // We found no match, or an ambiguous match and ran out of context.
            None
        }
        
    }
    
}
