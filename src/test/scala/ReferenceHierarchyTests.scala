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
        val pattern = "ACTAGT"
        val mappings: Seq[Option[Side]] = hierarchy.map(0, pattern)
        
        // Nothing should map because it's an ambiguous palindrome.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 0)
    }
    
    test("can map on level 1") {
        val pattern = "ACTAGT"
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern)
        
        mappings.foreach(println _)
        
        // All characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 6)
    }
    
    test("left-mapping answers are correct") {
        val pattern = "ACTAGT"
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern, Face.LEFT)
        
        println("Left mappings:")
        mappings.foreach(println _)
        
    }
    
    test("right-mapping answers are correct") {
        val pattern = "AATCTACTCC"
        val mappings: Seq[Option[Side]] = hierarchy.map(1, pattern, Face.RIGHT)
        
        println("Right mappings:")
        mappings.foreach(println _)
        
    }
    
}
