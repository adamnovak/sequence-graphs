import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.graph
import org.apache.spark.graph._
import edu.ucsc.soe.sequencegraph
import edu.ucsc.soe.sequencegraph._

//import org.broadinstitute.variant.vcf.VCFFileReader;

import org.rogach.scallop._;

object SequenceGraphs {
    def main(args: Array[String]) {
        
        // Option parsing code more or less stolen from various Scallop examples.
        val opts = new ScallopConf(args) {
            banner("""Usage: importVCF [OPTION]... vcfFile
                |Import a VCF file to sequence graph format.
                |Options:
                |""".stripMargin)
        
            val vcfFilename = trailArg[String]("vcfFilename", required = true)
            
            val version = opt[Boolean]("version", noshort = true, 
                descr = "Print version")
            val help = opt[Boolean]("help", noshort = true, 
                descr = "Show this message")

        } 
        
        println(opts.vcfFilename)
        
        /*    
        val vcfFile = File(args
    
        final VCFFileReader fileReader = new VCFFileReader(file);
        final VCFHeader fileHeader = fileReader.getFileHeader();
        */

    
    }
}

