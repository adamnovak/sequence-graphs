import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.graph
import org.apache.spark.graph._
import edu.ucsc.soe.sequencegraph
import edu.ucsc.soe.sequencegraph._

import java.io.File
import java.nio.file.Paths

import edu.unc.genomics.io.VCFFileReader

import org.rogach.scallop._;

object SequenceGraphs {
    def main(args: Array[String]) {
        
        // Option parsing code more or less stolen from various Scallop
        // examples.
        val opts = new ScallopConf(args) {
            guessOptionName = true
            banner("""Usage: importVCF [OPTION]... vcfFile
                |Import a VCF file to sequence graph format.
                |Options:
                |""".stripMargin)
            
            val vcfFile = trailArg[String](required = true)
            
            val version = opt[Boolean](noshort = true, 
                descr = "Print version")
            val help = opt[Boolean](noshort = true, 
                descr = "Show this message")

        } 
        
        // Get the actual String or die trying, and then make a Path from it,
        // and then open that.
        val vcfFile = new VCFFileReader(Paths.get(opts.vcfFile.get.get))
        
        println(vcfFile)
        
        
        /*    
        val vcfFile = File(args
    
        final VCFFileReader fileReader = new VCFFileReader(file);
        final VCFHeader fileHeader = fileReader.getFileHeader();
        */

    
    }
}

