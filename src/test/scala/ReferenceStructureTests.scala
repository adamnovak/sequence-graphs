package edu.ucsc.genome

import org.scalatest._

/**
 * Tests for making and using ReferenceStructure objects.
 */
class ReferenceStructureTests extends FMDSuite {
    
    // We will keep a StringReferenceStructure around
    var stringReference: StringReferenceStructure = null

    test("StringReferenceStructure can be created") {
        stringReference = new StringReferenceStructure(basename)
    }
    
    test("maps all bases in contig") {
    
        val pattern = "AATCTACTGC"
        
        val leftMappings = stringReference.mapFace(pattern, Face.LEFT)
        val rightMappings = stringReference.mapFace(pattern, Face.RIGHT)
        
        // Zip them together and disambiguate each pair. Note that (a, b).zipped
        // is of a type that provides a map that takes binary functions, while
        // a.zip(b).map takes only unary functions.
        val mappings = (leftMappings, rightMappings).zipped
            .map(stringReference.disambiguate(_, _))
        
        // All 10 characters ought to map.
        val mapped = mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum
        assert(mapped === 10)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the left, since map has
                // left-mapping semantics for its output.
                assert(mapping.face === Face.LEFT)
                // And on this first contig
                assert(mapping.coordinate < stringReference.getIndex
                    .getContigLength(0))
            }
            case None => {}
        }
        
    }
    
    test("maps all bases in reverse complement") {
    
        // This is the second contig.
        val pattern = "GCTAGTAGCTT"
        val mappings: Seq[Option[Side]] = stringReference.map(pattern)
        
        // The first 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 11)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the right, since map
                // has left-mapping semantics for its output and we're mapping
                // to a reverse strand.
                assert(mapping.face === Face.RIGHT)
                // And on this second contig
                assert(mapping.coordinate >= stringReference.getIndex
                    .getContigLength(0))
                assert(mapping.coordinate < stringReference.getIndex
                    .getContigLength(0) + stringReference.getIndex
                    .getContigLength(1))
            }
            case None => {}
        }
        
    }
    
    test("discards side-ambiguous mappings") {
        // Splice a couple sequences together on that T at base 5
        val pattern = "AATCTAGTAGCTT"
        
        val mapping = stringReference.map(pattern)(4)
        assert(mapping == None)
    }
}
