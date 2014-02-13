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
     * Test if the given Position is within this range.
     */
    def contains(position: Position): Boolean = {
        // It's in range if it's on the correct contig, and it's between start
        // and end inclusive. Face doesn't matter since we implicitly represent
        // the outer faces of the bases we end at.
        position.contig == contig && position.base >= start && 
            position.base <= end
    }
    
    /**
     * Test if this range is completely between the given positions, and does
     * not overlap them.
     */
    def between(position1: Position, position2: Position): Boolean = {
        if(position1.contig != contig || position2.contig != contig) {
            // We can't be between those positions, they're on some other
            // contig.
            false
        } else if (position1.base > position2.base) {
            // Flip them around
            between(position2, position1)
        } else if(position1.base < start && position2.base > end) {
            // We actually are between these positions
            true
        } else {
            // We aren't between these positions
            false
        }
    
    }
    
    /**
     * Output in contig:start-end format.
     */
    override def toString = "%s:%d-%d".format(contig, start, end)    
    
}
