package edu.ucsc.genome

import org.scalatest._

/**
 * Tests for making and using ReferenceHierarchy objects.
 */
class ReferenceHierarchyTests extends HierarchySuite {
    
    // We will keep a ReferenceHierarchy around
    var hierarchy: ReferenceHierarchy = null

    test("can be created") {
        // Load the hierarchy from the filename HierarchySuite feeds us.
        hierarchy = new ReferenceHierarchy(indexName)
    }
    
    test("can map on level 0") {
        // This is just seq1
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Side]] = hierarchy.map(0, pattern)
        
        // All characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 10)
    }
    
    test("maps merged bases to the same places") {
        // Have our two sequences.
        val seq1 = "AATCTACTGC"
        val seq2 = "AAGCTACTAGC"
        
        // The "C" at the start of "CTAC" ought to map to the same place in
        // both.
        val mapping1 = hierarchy.map(1, seq1)(3)
        val mapping2 = hierarchy.map(1, seq2)(3)
        
        assert(mapping1 != None)
        assert(mapping1 === mapping2)
        
    }
    
    test("can map on level 1") {
        // This is just seq1
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern)
        
        mappings.foreach(println _)
        
        // All characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 10)
    }
    
}
