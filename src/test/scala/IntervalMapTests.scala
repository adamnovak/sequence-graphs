package edu.ucsc.genome

import org.scalatest._
import fi.helsinki.cs.rlcsa.RangeVectorIterator

/**
 * Tests for the IntervalMap class, to make sure it doesn't screw up containment
 * or something.
 */
class IntervalMapTests extends FunSuite {

    var map: IntervalMap[String] = null

    test("can be created") {
        map = new IntervalMap[String]
    }
    
    test("can insert and merge") {
        map(9, 10) = "Hello"
        map(5, 6) = "Hello"
        map(7, 8) = "Hello"
        map(1, 4) = "World"
    }
    
    test("can retrieve exactly by merged interval") {
        assert(map.get(5, 10) === Some("Hello"))
        assert(map.get(1, 4) === Some("World"))
        
        assert(map(5, 10) === "Hello")
        assert(map(1, 4) === "World")
    }
    
    test("can retrieve by subinterval") {
        assert(map.get(6, 6) === Some("Hello"))
        assert(map.get(2, 4) === Some("World"))
        
        assert(map(7, 8) === "Hello")
        assert(map(1, 2) === "World")
    }
    
    test("ignores non-subintervals") {
        assert(map.get(0, 10) === None)
        assert(map.get(4, 6) === None)
        assert(map.get(-100, 100) === None)
    }
    
    test("ignores out-of-range intervals") {
        assert(map.get(0, 0) === None)
        assert(map.get(100, 100) === None)
    }
    
    test("correctly creates the value array") {
        val array = map.valueArray
        
        assert(array(0) === None)
        assert(array(1) === Some("World"))
        assert(array(2) === Some("Hello"))
        assert(array.size === 3)
    }
    
    test("correctly creates the range vector") {
        val iterator = new RangeVectorIterator(map.rangeVector)
        
        assert(iterator.rank(0) === 0)
        assert(iterator.rank(1) === 1)
        assert(iterator.rank(4) === 1)
        assert(iterator.rank(5) === 2)
        assert(iterator.rank(10) === 2)
        // This one is past the end
        assert(iterator.rank(11) === 3)
        // Select at the end here is strange (select(3) = 5), but we'll never
        // use it, so we don't bother testing it.
        
    }
    
    test("works when empty") {
        val emptyMap = new IntervalMap[String]
        assert(emptyMap.get(1, 3) === None)
        assert(emptyMap.valueArray.size === 0)
        
        val iterator = new RangeVectorIterator(emptyMap.rangeVector)
        // TODO: this probably ought to properly be 1, as 0 is past the end of a
        // truly empty vector, but we can't make a RangeVector that doesn't have
        // at least one space (which is implicitly a 0) or the vanilla RLCSA
        // code crashes.
        assert(iterator.rank(0) === 0)
        // This is past the end, and thus has rank 1 on an empty vector.
        assert(iterator.rank(100) === 1)
    } 


}
