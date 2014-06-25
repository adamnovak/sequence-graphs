package edu.ucsc.genome

import org.scalatest._

/**
 * Tests for making and using ReferenceHierarchy objects.
 */
class ReferenceHierarchyTests extends HierarchySuite {
    
    // We will keep a ReferenceHierarchy around
    var hierarchy: ReferenceHierarchy = null

    // Demand a palindrome.
    override def sequences = Seq("GATTACA", "GATTACA")

    test("can be created") {
        // Load the hierarchy from the filename HierarchySuite feeds us.
        hierarchy = new ReferenceHierarchy(indexName)
    }
    
    test("can't map on level 0 due to duplication") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.map(0, pattern)
        
        // Nothing should map because it's in there twice.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 0)
    }
    
    test("can left map on level 1 after merge") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern, Face.LEFT)
        
        // All characters ought to map, except the leftmost.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 6)
        
        assert(mappings(0) === None)
    }
    
    test("can right map on level 1 after merge") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern, Face.RIGHT)
        
        // All characters ought to map, except the rightmost.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 6)
        
        assert(mappings(6) === None)
    }
    
    test("can map on level 1 after merge") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern)
        
        // All characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 7)
    }
    
}
