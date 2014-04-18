package edu.ucsc.genome

import scala.util.parsing.combinator.RegexParsers

/**
 * Represents an interval in a linear contig. Uses 1-based inclusive indexing
 * and parses/outputs contig:start-end notation. Semantically, runs from the
 * left Side of the start base to the right Side of the end base. Start must be
 * less than or equal to end. Contig may not contain ":".
 *
 * TODO: Un-serializable this. We only have it serializable because Kryo argues
 * with scalatest and refuses to serialize it in tests.
 */
class BaseRange(val contig: String, val start: Long, val end: Long) 
    extends Serializable {
    
    if(start > end) {
        // We can't have backwards bounds
        throw new Exception("Can't create a negative-length range %s:%d-%d"
            .format(contig, start, end))
    }
    
    if(contig.contains(":")) {
        // We can't have pathological contig names.
        throw new Exception(
            "Cannot make a range on %s because contig name contains \":\""
            .format(contig))
    }
    
    def this(string: String) {
        // If only I could define some vals here this would be much simpler. And
        // we could have error checking. Everything would be great.
        this(string.split(":")(0), string.split(":")(1).split("-")(0).toLong, 
            string.split(":")(1).split("-")(1).toLong)
    }
    
    /**
     * Output in contig:start-end format.
     */
    override def toString = "%s:%d-%d".format(contig, start, end)    
    
}
