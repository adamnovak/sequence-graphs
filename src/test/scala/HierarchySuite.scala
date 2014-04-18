package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._
import scala.sys.process._

/**
 * A test utility class that provides a createIndex-generated hierarchy fixture
 * on disk.
 */
abstract class HierarchySuite extends FunSuite with BeforeAndAfterAll {

    // Holds the temporary directory
    var scratch: Path = null
    
    // Holds the directory of the index
    var indexName: String = null
    
    override def beforeAll = {
        // Make a new temp directory
        scratch = Files.createTempDirectory("test")
        
        // Decide on the index directory
        indexName = scratch.resolve("index").toString
    
        // Pick a name
        val fasta = scratch.resolve("fasta1").toString
        
        // Write a FASTA
        val fastaWriter = new FileWriter(fasta)
        
        // Write a sequence (Benedict's example)
        fastaWriter.write(">seq1\n")
        fastaWriter.write("AATCTACTGC\n")
        
        fastaWriter.close
        
        // Pick a name
        val fasta2 = scratch.resolve("fasta2").toString
        
        // Write a FASTA
        val fastaWriter2 = new FileWriter(fasta2)
        
        // Write a sequence (Benedict's example)
        fastaWriter2.write(">seq2\n")
        fastaWriter2.write("AAGCTACTAGC\n")
        
        fastaWriter2.close
        
        // Invoke the createIndex tool on the two FASTAs.
        Seq("./createIndex.sh", indexName, fasta, fasta2)!
        
    }
    
    override def afterAll = {
        // Clean it up
        FileUtils.deleteDirectory(new File(scratch.toString))
    }
    
}
