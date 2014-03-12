package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._

/**
 * Tests for making and using FMDIndex objects.
 */
class FMDIndexTests extends RLCSASuite {

    var index: FMDIndex = null;

    test("FMDIndex can be created") {
        index = new FMDIndex(basename)
    }
    
    test("getPosition works at start of forward strand") {
        val position = index.getPosition(0, 0, 1)
        
        assert(position.contig == "seq1")
        assert(position.base == 1)
        assert(position.face == Face.LEFT)
    }
    
    test("getPosition works at start of a reverse strand") {
        val position = index.getPosition(1, 0, 1)
        assert(position.contig == "seq1")
        assert(position.base == 10)
        assert(position.face == Face.RIGHT)
    }
    
    test("getPosition works at end of a forward strand") {
        val position = index.getPosition(0, 9, 1)
        
        assert(position.contig == "seq1")
        assert(position.base == 10)
        assert(position.face == Face.LEFT)
    }
    
    test("getPosition works at end of a reverse strand") {
        val position = index.getPosition(1, 9, 1)
        
        assert(position.contig == "seq1")
        assert(position.base == 1)
        assert(position.face == Face.RIGHT)
    }
    
    test("getPosition works on a whole strand") {
        val position = index.getPosition(0, 0, 10)
        
        assert(position.contig == "seq1")
        assert(position.base == 10)
        assert(position.face == Face.LEFT)
    }
    
    test("does not find missing substring") {
        assert(index.count("AAAAA") == 0)
    }
    
    test("finds entire contig") {
        assert(index.count("AATCTACTGC") == 1)
    }
    
    test("finds reverse complement") {
        assert(index.count("GCTAGTAGCTT") == 1)
    }
    
    test("finds common subsequence") {
        assert(index.count("AA") == 2)
    }
    
    test("fails to left-map left-ambiguous thing") {
        val mapping = index.leftMap("AATCTACTGC", 2)
        assert(mapping == None)
    }
    
    test("right-maps left-ambiguous thing") {
        val position = index.rightMap("AATCTACTGC", 2).get
        
        assert(position.contig == "seq1")
        assert(position.base == 2)
        assert(position.face == Face.RIGHT)
    }
    
    test("left-maps last base of contig") {
        val position = index.leftMap("AATCTACTGC", 10).get
        
        assert(position.contig == "seq1")
        assert(position.base == 10)
        assert(position.face == Face.LEFT)
    }
    
    test("left-maps last base of reverse complement") {
        val position = index.leftMap("GCTAGTAGCTT", 11).get
        
        assert(position.contig == "seq2")
        assert(position.base == 1)
        assert(position.face == Face.RIGHT)
    }
    
    test("left-maps unambiguous bases in contig") {
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Position]] = index.leftMap(pattern)
        
        // The last 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum == 8)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the left
                assert(mapping.face == Face.LEFT)
                // And on this contig
                assert(mapping.contig == "seq1")
            }
            case None => {}
        }
        
        
    }
    
    test("right-maps unambiguous bases in contig") {
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Position]] = index.rightMap(pattern)
        
        // The first 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum == 8)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the right
                assert(mapping.face == Face.RIGHT)
                // And on this contig
                assert(mapping.contig == "seq1")
            }
            case None => {}
        }
        
        
    }
    
    test("maps all bases in contig") {
    
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Position]] = index.map(pattern)
        
        // The first 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum == 10)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the left, since map has
                // left-mapping semantics for its output.
                assert(mapping.face == Face.LEFT)
                // And on this contig
                assert(mapping.contig == "seq1")
            }
            case None => {}
        }
        
    }
    
    test("maps all bases in reverse complement") {
    
        val pattern = "GCTAGTAGCTT"
        val mappings: Seq[Option[Position]] = index.map(pattern)
        
        // The first 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum == 11)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the right, since map
                // has left-mapping semantics for its output and we're mapping
                // to a reverse strand.
                assert(mapping.face == Face.RIGHT)
                // And on this contig
                assert(mapping.contig == "seq2")
            }
            case None => {}
        }
        
    }
    
    test("discards side-ambiguous mappings") {
        // Splice a couple sequences together on that T at base 5
        val pattern = "AATCTAGTAGCTT"
        
        val mapping = index.map(pattern, 5)
        assert(mapping == None)
    }
    
}
