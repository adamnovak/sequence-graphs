package edu.ucsc.genome

import org.scalatest._

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
    
    test("has sane ranges") {
        // Get all the BWT rows.
        val bwt = reference.getIndex.bwtTable
        // Organize into tuples of a Side and all its BWT rows.
        val sideContexts = reference.getRanges.map {
            case ((first, last), side) => 
            (side, (first until last + 1).map {
                (bwtIndex: Long) => bwt(bwtIndex.toInt)
            })
        }
        
        sideContexts.foreach { case (side, contexts) =>
            // Make sure they all start with the same letter.
            assert(contexts.map(_(0)).distinct.size === 1)
        }
    }
    
    
}
