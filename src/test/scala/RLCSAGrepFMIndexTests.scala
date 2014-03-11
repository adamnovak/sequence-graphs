package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._

/**
 * Tests for making and using RLCSAGrepFMIndex objects.
 */
class RLCSAGrepFMIndexTests extends RLCSASuite {

    var index: RLCSAGrepFMIndex = null;

    test("RLCAGrepFMIndex can be created") {
        index = new RLCSAGrepFMIndex(basename)
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
    
    test("fails to map ambiguous thing") {
        val mapping = index.map("AATCTACTGC", 2)
        assert(mapping == None)
    }
    
    test("maps entire contig") {
        val position = index.map("AATCTACTGC", 10).get
        
        assert(position.contig == "seq1")
        assert(position.base == 10)
        assert(position.face == Face.LEFT)
    }
    
    test("maps reverse complement") {
        val position = index.map("GCTAGTAGCTT", 11).get
        
        assert(position.contig == "seq2")
        assert(position.base == 1)
        assert(position.face == Face.RIGHT)
    }
    
}
