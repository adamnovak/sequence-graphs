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
    
    test("can map on level 1") {
        // This is just seq1
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern)
        
        // All characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 10)
    }
    
}
