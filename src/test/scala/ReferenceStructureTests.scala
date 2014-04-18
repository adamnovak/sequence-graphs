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
    
        println("====Starting Problematic Test====")
    
        val pattern = "AATCTACTGC"
        
        val leftMappings = stringReference.getIndex.map(pattern, Face.LEFT)
        val rightMappings = stringReference.getIndex.map(pattern, Face.RIGHT)
    
        // Zip them together and disambiguate each pair. Note that (a, b).zipped
        // is of a type that provides a map that takes binary functions, while
        // a.zip(b).map takes only unary functions.
        val mappings = (leftMappings, rightMappings).zipped  map(stringReference.disambiguate(_, _))
        
        //val mappings: Seq[Option[Position]] = stringReference.map(pattern)
        
        // All 10 characters ought to map.
        val mapped = mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum
        
        if(mapped != 10) {
            // Here is where we have problems.
            println("!!! Problem detected")
            
            // Wait for the user to acknowledge this
            System.in.read
            
            println("!!! Problem detected")
            
            println("Incorrect mappings:")
            println(mappings.mkString("\n"))
            
            println("Left was:")
            println(leftMappings.mkString("\n"))
            
            println("Right was:")
            println(rightMappings.mkString("\n"))
            
            // Try again
            val mappings2 = stringReference.map(pattern)
            
            val mapped2 = mappings2.map {
                case Some(_) => 1
                case None => 0
            }.sum
            
            println("On second attempt: %d/10 map:".format(mapped2))
            
            println(mappings2.mkString("\n"))
            
            // Wait for the user to acknowledge this
            System.in.read
        }
        
        assert(mapped === 10)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the left, since map has
                // left-mapping semantics for its output.
                assert(mapping.face === Face.LEFT)
                // And on this contig
                assert(stringReference.getIndex.positionToContigNumber(
                    mapping.coordinate) === 0)
            }
            case None => {}
        }
        
    }
    
    test("maps all bases in reverse complement") {
    
        val pattern = "GCTAGTAGCTT"
        val mappings: Seq[Option[Side]] = stringReference.map(pattern)
        
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
                assert(stringReference.getIndex.positionToContigNumber(
                    mapping.coordinate) === 1)
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
        collapsedReference = new CollapsedReferenceStructure(stringReference)
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
    
    test("CollapsedReferenceStructure maps previously ambiguous things") {
        // This would be ambiguous if it weren't for the merging.
        val mappings = collapsedReference.map("AA")
        
        println(mappings.mkString("\n"))
        
        assert(mappings(0) != None)
    }
    
}
