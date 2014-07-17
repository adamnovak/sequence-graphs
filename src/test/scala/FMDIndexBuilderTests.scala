package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._

import org.ga4gh.FMDIndexBuilder

/**
 * Tests for making and using FMDIndexBuilder objects.
 */
class FMDIndexBuilderTests extends FunSuite with BeforeAndAfterAll {

    // Holds the temporary directory
    var scratch: Path = null
    
    // Hold the object under test
    var builder: FMDIndexBuilder = null
    
    override def beforeAll = {
        // Make a new temp directory
        scratch = Files.createTempDirectory("test")
    }
    
    override def afterAll = {
        // Clean it up
        FileUtils.deleteDirectory(new File(scratch.toString))
    }
    
    // Run the tests
    
    test("FMDIndexBuilder can be created") {
        // Decide on the index filename
        val basename = scratch.resolve("index").toString
    
        // Make the builder
        builder = new FMDIndexBuilder(basename)
    }
    
    test("FMDIndexBuilder can add first sequence") {
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
    
    test("FMDIndexBuilder can add additional sequences") {
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
    
    test("FMDIndexBuilder can finish index") {
        builder.build
    }
    
    // TODO: Somehow check this index without relying on the index readers to be
    // correct.
}
