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
                
            val directory = trailArg[String](required = true,
                descr = "Directory to save under (will overwrite)")
            
            val dotFile = opt[String](
                descr = "Save a GraphViz graph to this .dot file")
            
            val version = opt[Boolean](noshort = true, 
                descr = "Print version")
            val help = opt[Boolean](noshort = true, 
                descr = "Show this message")

        } 
        
        // What sample are we importing?
        val sample = opts.sampleName.get.get
        
        // Make a new SequenceGraphBuilder to build its graph.
        val graph = new SequenceGraphBuilder(sample, "reference")
        
        // Get the actual File or die trying, and then make VcfParser parse
        // that. We need to invoke the VcfParser functor first for some reason.
        VcfParser().parseFile(opts.vcfFile.get.get, false) { 
            (info, entries) =>
            // We get the VCF metadata and an iterator of VCF entries.
            // Import our sample
            importSample(graph, info, entries, sample)
        }
        
        opts.dotFile.get map { (file) =>
            // Write out a graphviz graph
            graph.writeDotFile(file)
        }
        
        // Write Parquet files to the current directory.
        graph.writeParquetFiles(opts.directory.get.get)
        
        // TODO: Flush logging output from Parquet/Spark stuff.
        println("VCF imported")
        
    }
    
    /**
     *
     *   Import a sample into the given SequenceGraphBuilder from a VCF, creating
     *   a list of Sides and a list of SequenceGraphEdges, as well as Sites and
     *   Breakpoints for them to be in, and Alleles they contain.
     *   
     *   Takes the VCF header metadata, an iterator over the VCF's variants, and the
     *   name of the sample to import.
     * 
     */
    def importSample(builder: SequenceGraphBuilder, info: VcfInfo, 
        entries: Iterator[(Variant, List[Metadata.Format], 
        List[List[List[VcfValue]]])], sampleName: String) {
        // TODO: Can I do something to not have to include this big ugly type?
        // The vcfimp flatten tool accomplishes this by returning an anonymous
        // function typed VcfParser.Reader[Either[String, Unit]] which lets the
        // types on the internal anonymous function be inferred, but I don't
        // like the wierd partial application format that gives me.
        
        
        // Get the index of the sample we want
        val sampleIndex = info.samples indexWhere { (sample) => 
            sample.id.id == sampleName
        }
        
        if(sampleIndex < 0) {
            throw new Exception("Sample %s not found in file".format(
                sampleName))
        } else {
            println("Importing sample %s at index %d...".format(
                sampleName, sampleIndex))
        }
        
        // Was the last variant processed phased?
        var lastCallPhased: Boolean = true
        
        // What contig was the last variant on? We assume the variants are
        // sorted by contig and position.
        var lastContig: String = null
        
        // Where did the last variant end on that contig?
        var lastEnd: Int = 0
        
        for((variant, formats, samples) <- entries) {
            // For each variant in the file
            
            // What contig is it on?
            val contig: String = variant.chromosome match {
                case Left(vcfid) => vcfid.toString
                case Right(string) => "chr" + string
            }
            
            // Where does it start in the reference?
            val referenceStart = variant.position
            
            // Where does it end in the reference?
            val referenceEnd = referenceStart + variant.reference.size
            
            // What alleles are available? The reference and all the alternates.
            val alleles = Right(variant.reference) :: variant.alternates
            
            if(contig != lastContig && lastContig != null) {
                // We're ending an old contig and starting a new one.
                
                // Add trailing telomeres
                builder.close(lastContig, 0)
                builder.close(lastContig, 1)
                
                // Create new leading telomeres
                builder.getLastSide(contig, 0)
                builder.getLastSide(contig, 1)
                
                // Record that the leading telomeres are phased, and on this
                // contig, and at the very start.
                lastCallPhased = true
                lastContig = contig
                lastEnd = 0
            }
            
            // This holds the geontype values read for the sample
            val sampleValues = samples(sampleIndex)
            
            // Find the index of the genotype field
            val genotypeIndex = formats indexWhere { (format) =>
                format.id.id == "GT"
            }
            
            // How far are we from the last variant? If this is 0, no Anchors
            // will actually get added.
            val referenceDistance = referenceStart - lastEnd
            
            if(genotypeIndex != -1) {
                // We actually have a GT field for this variant
            
                // Pull out the genotype string as a String.
                val genotypeString = sampleValues(genotypeIndex).head match {
                    case VcfString(value) => value
                    case _ => throw new Exception("Got a non-string GT value")
                }
            
                // At the moment we support only phased or unphased diploid
                // genomes.
                
                // Is the string a phased genotype?
                val phased = genotypeString contains "|"
                
                // Add the intermediate Anchors between the last variant site
                // and this one
                if(phased && lastCallPhased) {
                    // We need phased anchors, since both this call and the
                    // previous one are phased
                    
                    builder.addAnchor(contig, List(0), referenceDistance)
                    builder.addAnchor(contig, List(1), referenceDistance)
                    
                } else {
                    // We need an unphased anchor.
                    
                    builder.addAnchor(contig, List(0, 1), referenceDistance)
                }
                
                // TODO: Add a backwards empty AlleleGroup to discard phasing if
                // there is no room for an Anchor to do it.
                
                // Now we want to add the actual AlleleGroups.
                
                // Pull out and intify the two allele indices
                val alleleIndices = genotypeString.split("[|/]").map(_.toInt)
                
                // TODO: Map over phases.
                
                // The alternates are in the variant.alternates list.
                // Load the allele in phase 0
                val allele0: Either[Either[Breakend, VcfId], String] = 
                    alleles(alleleIndices(0))
                    
                // Load the allele in phase 1
                val allele1: Either[Either[Breakend, VcfId], String] = 
                    alleles(alleleIndices(1))
                
                // For now we handle only string alleles.
                
                val alleleString0 = allele0 match {
                    case Right(string) => string
                    case Left(_) => throw new Exception(
                        "Breakends and ID'd alleles not supported")
                }
                
                val alleleString1 = allele1 match {
                    case Right(string) => string
                    case Left(_) => throw new Exception(
                        "Breakends and ID'd alleles not supported")
                }
                
                // Add the Alleles to the appropriate Phases. TODO: if they are
                // the same and unphased, add just one AlleleGroup
                builder.addAllele(contig, List(0), new Allele(alleleString0), 
                    referenceEnd - referenceStart)
                builder.addAllele(contig, List(1), new Allele(alleleString1), 
                    referenceEnd - referenceStart)
                    
                // Save information about this site
                lastCallPhased = phased
                lastContig = contig
                lastEnd = referenceEnd
            }
            
        }
        
        
    }
}





















