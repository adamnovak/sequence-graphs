package edu.ucsc.genome

import org.scalatest._
import org.ga4gh.FMDIndex

/**
 * Tests for making and using MergedReferenceStructure objects.
 */
class MergedReferenceStructureTests extends HierarchySuite {
    
    // Have a merged reference structure
    var reference: MergedReferenceStructure = null

    test("MergedReferenceStructure can be created") {
        reference = new MergedReferenceStructure(new FMDIndex(indexName + 
            "/index.basename"), indexName + "/level1")
    }
    
}
