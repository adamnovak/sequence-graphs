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
        // Dump all the bits and the BWT rows they go with.
        (reference.getIndex.bwtTable zip reference.bits).zipWithIndex.map(println _)
        
        // Dump all the letters
        (reference.getIndex.firstColumn zip reference.getIndex.lastColumn zip reference.getIndex.suffixes).zipWithIndex.map(println _)
        
        // Dump all the ranges in order
        reference.getRanges.map(println _)
        
        // And collected by Side
        reference.sideArray.distinct.map { side =>
            println("%s: %s".format(side, reference.getRanges(side)
                .mkString(", ")))
        }
    }
    
    
}
