package edu.ucsc.genome

import org.scalatest._

/**
 * Tests for making and using ReferenceHierarchy objects on a plaindrome.
 */
class ReferenceHierarchyPalindromeTests extends HierarchySuite {
    
    // We will keep a ReferenceHierarchy around
    var hierarchy: ReferenceHierarchy = null

    // Demand a palindrome.
    override def sequences = Seq("ACTAGT")

    test("can be created") {
        // Load the hierarchy from the filename HierarchySuite feeds us.
        hierarchy = new ReferenceHierarchy(indexName)
    }
    
    test("can't map on level 0 due to palindrome") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.mapLevel(0, pattern)
        
        // Nothing should map because it's an ambiguous palindrome.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 0)
    }
    
    test("can't map on level 1 either because merging is now map-based") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.mapLevel(1, pattern)
        
        // All characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 0)
    }
    
    test("left-mapping is reverse of right-mapping") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.mapFaceLevel(1, pattern, 
            Face.LEFT)
        
        val mappings2: Seq[Option[Side]] = hierarchy.mapFaceLevel(1, pattern, 
            Face.RIGHT)
        
        assert(mappings === mappings2.reverse)
        
    }
}
