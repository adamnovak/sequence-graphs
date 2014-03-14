package edu.ucsc.genome

import java.util.{TreeMap, NoSuchElementException}
import fi.helsinki.cs.rlcsa._
import scala.collection.JavaConversions._
import scala.collection.mutable.ArrayBuffer

/**
 * Represents a mutable map from non-overlapping intervals of longs to some
 * other type. Given a query interval, can get the entry corresponding to the
 * interval that completely contains the query, or None if no interval
 * completely contains the query.
 *
 * Intervals are inclusive.
 *
 * Used for mapping to a reference structure using an index of a lower-level
 * (phased) graph. Could be implemented efficiently as a binary tree or as a bit
 * vector marking interval ends and an array of result values.
 */
class IntervalMap[ValueType] extends Serializable {
    
    
    // We keep our intervals in this sorted NavigableMap from interval start to
    // (interval end, value). We then can get the interval starting before the
    // start of a query, and just check that it covers the end.
    val map = new TreeMap[Long, (Long, ValueType)]
    
    /**
     * Read the value for the interval covering this one, or None if no interval
     * covers this one.
     */
    def get(interval: (Long, Long)): Option[ValueType] = {
        
        // Get the highest interval start <= the start of the query
        var entry = map.floorEntry(interval._1)
        
        if(entry != null) {
            // Unpack the interval we found
            val (end, value) = entry.getValue
            
            if(end >= interval._2) {
                // The interval we found covers our entire query. We have a
                // match.
                Some(value)
            } else {
                // We found an itnerval but the query fell off the end.
                None
            }
        } else {
            // We didn't find any interval starting before the query.
            None
        }
    }
    
    /**
     * Syntactic sugar for get without doubled parentheses.
     */
    def get(start: Long, end: Long): Option[ValueType] = {
        get((start, end))
    }
    
    /**
     * Get the value for the interval covering the viven interval. Throws an
     * exception if no such value exists.
     */
    def apply(interval: (Long, Long)): ValueType = {
        get(interval) match {
            // Return it if you find it
            case Some(thing) => thing
            // Otherwise complain
            case None => throw new NoSuchElementException(
                "key not found: %s".format(interval))
        }
    }
    
    /**
     * Syntactic sugar for apply without doubled parentheses.
     */
    def apply(start: Long, end: Long): ValueType = {
        apply((start, end))
    }
    
    /**
     * Store the given value at the given interval. It must not overlap with any
     * other interval.
     */
    def update(interval: (Long, Long), value: ValueType): Unit = {
        if(interval._2 < interval._1) {
            // Ban backwards intervals
            throw new Exception("Cannot insert backwards interval %d to %d"
                .format(interval._1, interval._2))
        }
        
        // Just assume we don't overlap, for speed. Stick the interval end and
        // the value in under the interval start.
        map.put(interval._1, (interval._2, value))
    }
    
    /**
     * Syntactic sugar for update without doubled parentheses.
     */
    def update(start: Long, end: Long, value: ValueType): Unit = {
        update((start, end), value)
    }
    
    /**
     * Turn into a RangeVector and sequence of values for ranges in that vector,
     * for efficient mapping.
     *
     * TODO: Cache the result.
     */
    def bake: (RangeVector, Seq[Option[ValueType]]) = {
        // Make a new encoder with this arbitrary block size.
        val encoder = new RangeEncoder(32)
        
        // Make an ArrayBuffer to store pointers to the things we map to by
        // RangeVector range number. Some "ranges" may really be the ranges
        // between ranges we have values for, so this needs to be full of
        // Options.
        val rangeMapping = new ArrayBuffer[Option[ValueType]]
        
        // Keep track of the end of the last interval, so you know if you
        // skipped anything.
        var lastEnd: Long = -1
        
        map.foreach { case (start, (end, value)) =>
            if(start != lastEnd + 1) {
                // We skipped something. Add a fake range at this number.
                rangeMapping += None
            }
            
            // Go through intervals (in order, since it's a sorted map).
            // Mark the start of each.
            encoder.addBit(start)
            // Mark the position after the end of each in case we skip
            // something.
            encoder.addBit(end + 1)
            
            // Say that this range belongs to this value
            rangeMapping += Some(value)
            
            // Save the endpoint of the interval
            lastEnd = end
        }
        
        // Finish and return the RangeVector and out range to value mapping.
        encoder.flush()
        (new RangeVector(encoder, lastEnd + 1), rangeMapping)
    }
    
    /**
     * Convert to a RangeVector for efficient mapping.
     */
    def rangeVector: RangeVector = bake._1
    
    /**
     * Get the array of values by interval index, for mapping back from ranges
     * in our RangeVector representation to the actual values in the map.
     */
    def valueArray: Seq[Option[ValueType]] = bake._2
}

