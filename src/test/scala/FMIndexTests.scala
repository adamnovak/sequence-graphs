package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._

/**
 * Tests for making and using FMIndex objects.
 */
class FMIndexTests extends FunSuite with BeforeAndAfterAll {

    // Holds the temporary directory
    var scratch: Path = null
    
    // These hold the objects under test
    var builder: RLCSABuilder = null
    var index: FancyFMIndex = null
    
    override def beforeAll = {
        // Make a new temp directory
        scratch = Files.createTempDirectory("test")
    }
    
    override def afterAll = {
        // Clean it up
        FileUtils.deleteDirectory(new File(scratch.toString))
    }
    
    // Run the tests
    
    test("RLCSABuilder can be created") {
        // Decide on the index filename
        val basename = scratch.resolve("index").toString
    
        // Make the builder
        builder = new RLCSABuilder(basename)
    }
    
    test("RLCSABuilder can add first sequence") {
        // Pick a name
        val fasta = scratch.resolve("fasta1").toString
        
        // Write a FASTA
        val fastaWriter = new FileWriter(fasta)
        
        // Write a sequence (Benedict's example)
        fastaWriter.write(">seq1\n")
        fastaWriter.write("AATCTACTGC\n")
        
        fastaWriter.close
        
        // Add it to the index
        builder.add(fasta)
    }
    
    test("RLCSABuilder can add additional sequences") {
        // Pick a name
        val fasta = scratch.resolve("fasta2").toString
        
        // Write a FASTA
        val fastaWriter = new FileWriter(fasta)
        
        // Write a sequence (Benedict's example)
        fastaWriter.write(">seq2\n")
        fastaWriter.write("AAGCTACTAGC\n")
        
        fastaWriter.close
        
        // Add it to the index
        builder.add(fasta)
    }
    
    test("RLCAGrepFMIndex can be created") {
        index = builder.getIndex
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
