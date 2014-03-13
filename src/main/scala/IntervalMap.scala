package edu.ucsc.genome

import java.util.{TreeMap, NoSuchElementException}

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
}

