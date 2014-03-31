package edu.ucsc.genome

import org.scalatest._

/**
 * Tests for making and using ReferenceHierarchy objects.
 */
class ReferenceHierarchyTests extends RLCSASuite with SparkSuite {
    
    // We will keep a ReferenceHierarchy around
    var hierarchy: ReferenceHierarchy = null

    test("ReferenceHierarchy can be created") {
        // Use the SparkContext sc from SparkSuite, and the index basename from
        // RLCSASuite.
        hierarchy = new ReferenceHierarchy(sc, new FMDIndex(basename), 
            Seq(Unmerged(), NonSymmetric(2)))
            
        // Build the graphs and hierarchy levels.
        hierarchy.initialize
    }
    
    test("has the appropriate number of levels") {
        // Should be the base contigs and the two levels we specified.
        assert(hierarchy.levels.size === 3)
    }
    
    
    
}
