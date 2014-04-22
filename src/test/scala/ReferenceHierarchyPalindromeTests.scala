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
        val mappings: Seq[Option[Side]] = hierarchy.map(0, pattern)
        
        // Nothing should map because it's an ambiguous palindrome.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 0)
    }
    
    test("can map on level 1") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern)
        
        mappings.foreach(println _)
        
        // All characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === sequences(0).size)
    }
    
    test("left-mapping is reverse of right-mapping") {
        val pattern = sequences(0)
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern, Face.LEFT)
        
        val mappings2: Seq[Option[Side]] = hierarchy.map(1, pattern, Face.RIGHT)
        
        println(mappings.mkString("\n"))
        println(mappings2.mkString("\n"))
        
        assert(mappings === mappings2.reverse)
        
    }
}
