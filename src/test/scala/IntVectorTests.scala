package edu.ucsc.genome

import org.ga4gh.IntVector
import org.scalatest._

/**
 * Tests for the IntVector, which is supposed to hold 64-bit signed integers. It
 * is really a C++ vector. Makes sure it works through the Java bindings.
 */
class IntVectorTests extends FunSuite {

    var vector: IntVector = null

    test("can be created") {
        vector = new IntVector
    }
    
    test("can add to it") {
        vector.add(-1)
        vector.add(2)
        vector.add(1000)
    }
    
    test("can read from it") {
        assert(vector.get(0) === -1)
        assert(vector.get(1) === 2)
        assert(vector.get(2) === 1000)
    }


}
