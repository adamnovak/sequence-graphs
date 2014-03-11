package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._

/**
 * A test utility class that provides an RLCSA FMD-index fixture on disk.
 */
abstract class RLCSASuite extends FunSuite with BeforeAndAfterAll {

    // Holds the temporary directory
    var scratch: Path = null
    
    // Holds the basename of the index
    var basename: String = null
    
    override def beforeAll = {
        // Make a new temp directory
        scratch = Files.createTempDirectory("test")
        
        // Decide on the index filename
        basename = scratch.resolve("index").toString
    
        // Make the builder
        val builder = new RLCSABuilder(basename)
        
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
        
        // Pick a name
        val fasta2 = scratch.resolve("fasta2").toString
        
        // Write a FASTA
        val fastaWriter2 = new FileWriter(fasta2)
        
        // Write a sequence (Benedict's example)
        fastaWriter2.write(">seq2\n")
        fastaWriter2.write("AAGCTACTAGC\n")
        
        fastaWriter2.close
        
        // Add it to the index
        builder.add(fasta2)
    }
    
    override def afterAll = {
        // Clean it up
        FileUtils.deleteDirectory(new File(scratch.toString))
    }
    
}
