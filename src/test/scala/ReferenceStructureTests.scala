package edu.ucsc.genome

import org.scalatest._

/**
 * Tests for making and using ReferenceStructure objects.
 */
class ReferenceStructureTests extends RLCSASuite {
    
    // We will keep a StringReferenceStructure around
    var stringReference: StringReferenceStructure = null
    // And a collapsed reference structure on top of that
    var collapsedReference: CollapsedReferenceStructure = null

    test("StringReferenceStructure can be created") {
        stringReference = new StringReferenceStructure(basename)
    }
    
    test("maps all bases in contig") {
    
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Position]] = stringReference.map(pattern)
        
        println("All bases in contig")
        println(mappings.mkString("\n"))
        
        // All 10 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 10)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the left, since map has
                // left-mapping semantics for its output.
                assert(mapping.face === Face.LEFT)
                // And on this contig
                assert(mapping.contig === "seq1")
            }
            case None => {}
        }
        
    }
    
    test("maps all bases in reverse complement") {
    
        val pattern = "GCTAGTAGCTT"
        val mappings: Seq[Option[Position]] = stringReference.map(pattern)
        
        // The first 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum === 11)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the right, since map
                // has left-mapping semantics for its output and we're mapping
                // to a reverse strand.
                assert(mapping.face === Face.RIGHT)
                // And on this contig
                assert(mapping.contig === "seq2")
            }
            case None => {}
        }
        
    }
    
    test("discards side-ambiguous mappings") {
        // Splice a couple sequences together on that T at base 5
        val pattern = "AATCTAGTAGCTT"
        
        val mapping = stringReference.map(pattern)(4)
        assert(mapping == None)
    }
    
    test("CollapsedReferenceStructure can be created on top") {
        collapsedReference = new CollapsedReferenceStructure(stringReference,
            "merged")
    }
    
    test("CollapsedReferenceStructure can have positions merged") {
        // Construct Benedict's example
        // Merge the initial run of identical bases
        collapsedReference.merge("seq1", 1, "seq2", 1)
        collapsedReference.merge("seq1", 2, "seq2", 2)
        // Pass each variant of the SNP
        collapsedReference.pass("seq1", 3)
        collapsedReference.pass("seq2", 3)
        // Merge the next 5 identical bases
        collapsedReference.merge("seq1", 4, "seq2", 4)
        collapsedReference.merge("seq1", 5, "seq2", 5)
        collapsedReference.merge("seq1", 6, "seq2", 6)
        collapsedReference.merge("seq1", 7, "seq2", 7)
        collapsedReference.merge("seq1", 8, "seq2", 8)
        // Pass the extra base in seq2
        collapsedReference.pass("seq2", 9)
        // Merge the last two identical bases (now at an offset)
        collapsedReference.merge("seq1", 9, "seq2", 10)
        collapsedReference.merge("seq1", 10, "seq2", 11)
    }
    
    test("CollapsedReferenceStructure maps to merged positions") {
        val mappings = collapsedReference.map("AATCTACTGC")
        
        println(mappings.mkString("\n"))
        
        mappings.foreach {
            case Some(position) => assert(position.contig === "merged")
            case None => Unit
        }
    }
    
    test("CollapsedReferenceStructure maps previously ambiguous things") {
        // This would be ambiguous if it weren't for the merging.
        val mappings = collapsedReference.map("AA")
        
        println(mappings.mkString("\n"))
        
        val position = mappings(0).get
        
        assert(position.contig === "merged")
        assert(position.base === 1)
        assert(position.face === Face.LEFT)
    }
    
    test("CollapsedReferenceStructure maps to passed positions") {
        // Seq1 has the T, so it gets base number 3.
        val position = collapsedReference.map("AATCT")(2).get
        
        assert(position.contig === "merged")
        assert(position.base === 3)
        assert(position.face === Face.LEFT)
    }
    
    
    
}
