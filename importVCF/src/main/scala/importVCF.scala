import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.graph
import org.apache.spark.graph._
import edu.ucsc.genome._

import java.io.File
import java.nio.file.Paths

// We want to read VCF files
import ca.innovativemedicine.vcf._
import ca.innovativemedicine.vcf.parsers._

// We want to be able to loop over Java iterators: load a bunch of conversions.
// See <http://stackoverflow.com/a/1625252/402891>
import scala.collection.JavaConversions._

// We want to parse command-line arguments
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
            
            // What file should we read?
            val vcfFile = trailArg[File](required = true,
                descr = "VCF file to open")
            // What sample should we import?
            val sampleName = trailArg[String](required = true,
                descr = "Sample to import")
            
            val version = opt[Boolean](noshort = true, 
                descr = "Print version")
            val help = opt[Boolean](noshort = true, 
                descr = "Show this message")

        } 
        
        
        // Get the actual File or die trying, and then make VcfParser parse
        // that. We need to invoke the VcfParser functor first for some reason.
        VcfParser().parseFile(opts.vcfFile.get.get, false) { 
            (metadata, entries) =>
            // We get the VCF metadata and an iterator of VCF entries.
            importSample(metadata, entries, opts.sampleName.get.get)
            
        }
        
    }
    
    /**
    
        Import a sample from a VCF, creating a list of Sides and a list of
        SequenceGraphEdges, as well as Sites and Breakpoints for them to be in,
        and Alleles they contain.
        
        Takes the VCF metadata, an iterator over the VCF's variants, and the
        name of the sample to import.
    
    */
    def importSample(metadata : VcfInfo, 
        entries : Iterator[(Variant, List[Metadata.Format], 
        List[List[List[VcfValue]]])], sampleName : String) {
        // TODO: Can I do something to not have to include this big ugly type?
        // The vcfimp flatten tool accomplishes this by returning an anonymous
        // function typed VcfParser.Reader[Either[String, Unit]] which lets the
        // types on the internal anonymous function be inferred, but I don't
        // like the wierd partial application format that gives me.
        
        // Make a SequenceGraphBuilder to build up our sample's graph
        val builder = new SequenceGraphBuilder(sampleName, "reference")
        
        println(metadata)
        println("NS is:")
        println(metadata.getTypedMetadata[Metadata.Info](VcfId("NS")))
        println("DP is:")
        println(metadata.getTypedMetadata[Metadata.Info](VcfId("DP")))
        println("AF is:")
        println(metadata.getTypedMetadata[Metadata.Info](VcfId("AF")))
        println("AA is:")
        println(metadata.getTypedMetadata[Metadata.Info](VcfId("AA")))
        println("HQ is:")
        println(metadata.getTypedMetadata[Metadata.Format](VcfId("HQ")))
        println("DP is:")
        println(metadata.getTypedMetadata[Metadata.Format](VcfId("DP")))
        
        for((variant, formats, sample) <- entries) {
            println(sample)
        }
        
        
    }
}





















