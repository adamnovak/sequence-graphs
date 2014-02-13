package edu.ucsc.genome

import org.scalatest._

/**
 * Tests for the BaseRange class, to make sure it doesn't screw up containment
 * or something.
 */
class BaseRangeTests extends FunSuite {

    var range: BaseRange = null

    test("can be parsed from a string") {
        range = new BaseRange("chr1:10-15")
        assert(range.contig == "chr1")
        assert(range.start == 10)
        assert(range.end == 15)
    }
    
    test("can be stringified") {
        assert(range.toString == "chr1:10-15")
    }

    test("contains things in the middle") {
        assert(range contains new Position("chr1", 12, Face.LEFT))
    }
    
    test("contains things on the left edge") {
        assert(range contains new Position("chr1", 10, Face.LEFT))
    }
    
    test("contains things on the right edge") {
        assert(range contains new Position("chr1", 15, Face.RIGHT))
    }
    
    test("doesn't contain things on other contigs") {
        assert(!(range contains new Position("chr2", 12, Face.LEFT)))
    }
    
    test("doesn't contain things off the left") {
        assert(!(range contains new Position("chr1", 9, Face.RIGHT)))
    }
    
    test("doesn't contain things off the right") {
        assert(!(range contains new Position("chr1", 16, Face.LEFT)))
    }
    
    test("is between what it should be between") {
        assert(range.between(new Position("chr1", 5, Face.RIGHT), 
            new Position("chr1", 20, Face.RIGHT)))
    }
    
    test("is not between things when the left one overlaps it") {
        assert(!range.between(new Position("chr1", 12, Face.RIGHT), 
            new Position("chr1", 20, Face.RIGHT)))
    }

}
