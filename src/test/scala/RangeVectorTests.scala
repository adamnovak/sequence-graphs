package edu.ucsc.genome

import fi.helsinki.cs.rlcsa._
import org.scalatest._

/**
 * Tests for the RangeVector, RangeVectorIterator, and RangeEncoder from RLCSA,
 * to make sure they work through the Java bindings.
 */
class RangeVectorTests extends FunSuite {

    var vector: RangeVector = null
    var encoder: RangeEncoder = null

    test("has constructable encoder") {
        // Make a new encoder with a block size of 32
        encoder = new RangeEncoder(32)
    }
    
    test("can add some runs") {
        // Runs are (start, length) and are 1s. You can skip as much as you
        // want, but must add them from left to right.
        encoder.addRun(1, 1)
        encoder.addRun(9, 2)
        encoder.addRun(19, 1)
        encoder.flush()
    }
    
    test("can build vector") {
        // Build with a "universe size" of 20, which is the total length of the
        // vector (in case there are 0s after the last 1).
        vector = new RangeVector(encoder, 20)
    }
    
    test("can get ranks") {
        val iterator = new RangeVectorIterator(vector)
        
        assert(iterator.rank(0) === 0)
        assert(iterator.rank(1) === 1)
        assert(iterator.rank(2) === 1)
        assert(iterator.rank(8) === 1)
        assert(iterator.rank(9) === 2)
        assert(iterator.rank(10) === 3)
        assert(iterator.rank(11) === 3)
        assert(iterator.rank(18) === 3)
        assert(iterator.rank(19) === 4)
    }
    
    test("rank past the end produces total count") {
        val iterator = new RangeVectorIterator(vector)
        
        assert(iterator.rank(100) === 4)
    }
    
    test("can select") {
        val iterator = new RangeVectorIterator(vector)
        
        assert(iterator.select(0) === 1)
        assert(iterator.select(1) === 9)
        assert(iterator.select(2) === 10)
        assert(iterator.select(3) === 19)
    }
    
    test("select out of range produces one-past-the-end") {
        val iterator = new RangeVectorIterator(vector)
        
        assert(iterator.select(500) === 20)
    }
    


}
