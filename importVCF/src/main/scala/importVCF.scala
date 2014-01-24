package edu.ucsc.genome.ImportVCF

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
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
import org.rogach.scallop._

// We want reasonable logging from Parquet.
import java.util.logging._
import parquet.Log

import org.apache.hadoop.mapreduce.Job

class DiscardHandler extends java.util.logging.Handler {
    def close() {}
    def flush() {}
    def publish(record: LogRecord) {}
}

object SequenceGraphs {
    def main(args: Array[String]) {
        
        // Option parsing code more or less stolen from various Scallop
        // examples.
        val opts = new ScallopConf(args) {
            guessOptionName = true
            banner("""Usage: importVCF [OPTION] {--dot-file dotFile | 
                 --parquet-dir parquetDir} vcfFile sample
                |Import a VCF file to sequence graph format.
                |Options:
                |""".stripMargin)
            
            // What file should we read?
            val vcfFile = trailArg[File](required = true,
                descr = "VCF file to open")
            // What sample should we import?
            val sampleName = trailArg[String](required = true,
                descr = "Sample to import")
            
            val chromSizes = opt[File](
                descr = "File of chromosome sizes to read, for end telomeres")
            
            val dotFile = opt[String](
                descr = "Save a GraphViz graph to this .dot file")
                
            val parquetDir = opt[String](
                descr = "Save Parquet files in this directory (will overwrite)")
            
            val version = opt[Boolean](noshort = true, 
                descr = "Print version")
            val help = opt[Boolean](noshort = true, 
                descr = "Show this message")

        } 
        
        // Parquet "helpfully" forcibly adds an INFO-level logger to console if
        // we don't configure its logger specifically (i.e. if we let its log
        // messages pass through to the root logger). It also overrides any
        // level we set on its logger, resetting it to INFO. See
        // <https://github.com/Parquet/parquet-mr/blob/master/parquet-
        // common/src/main/java/parquet/Log.java>
        
        // The solution is to add a warning-level console logger specifically
        // for Parquet, and tell it not to propagate messages up to the root
        // logger.
        
        // This could be done through the properties-file-based Java logging
        // configuration mechanism, but that would require telling Java how to
        // *find* the properties file, which in turn requires either setting a
        // JVM command-line option to a filesystem path or messing about with
        // the Java preferences API and/or resource streams. See
        // <http://stackoverflow.com/q/805701/402891>
        
        // Get the logger Parquet is going to use (before Parquet's static
        // logger initialization code can run)
        val parquetLogger = Logger.getLogger(
            classOf[Log].getPackage().getName())
        // Attach our own ConsoleHandler for it
        val consoleHandler = new ConsoleHandler()
        consoleHandler.setLevel(Level.WARNING)
        parquetLogger.addHandler(consoleHandler)
        // Tell it not to send messages up
        parquetLogger.setUseParentHandlers(false)
        
        // Parse the chromosome sizes file (<name>\t<int> TSV) into a map from
        // chromosome names to lengths
        val chromSizes = opts.chromSizes.get.map { (file) =>
            import scala.io.Source
            Source.fromFile(file).getLines.map { (line) =>
                // Split each line on the tab
                val parts = line.split("\t")
                // Make the second item an int
                (parts(0), parts(1).toInt)
            }.toMap
        }
        
        // What sample are we importing?
        val sample = opts.sampleName.get.get
        
        // Make a new SequenceGraphBuilder to build its graph, of the right type
        // to write whatever the user is asking for.
        val graph = if(opts.parquetDir.get isDefined) {
            // The user wants to write Parquet. TODO: what if they also asked
            // for a dot file?
            println("Writing to Parquet")
            
            new ParquetSequenceGraphBuilder(sample, "reference", 
                opts.parquetDir.get.get)
        } else if (opts.dotFile.get isDefined) {
            // The user wants to write GraphViz
            println("Writing to Graphviz")
            
            new GraphvizSequenceGraphBuilder(sample, "reference", 
                opts.dotFile.get.get)
        } else {
            // The user doesn't know what they're doing
            throw new Exception("Must specify --dot-file or --parquet-dir")
        }
        
        // Get the actual File or die trying, and then make VcfParser parse
        // that. We need to invoke the VcfParser functor first for some reason.
        VcfParser().parseFile(opts.vcfFile.get.get, false) { 
            (info, entries) =>
            // We get the VCF metadata and an iterator of VCF entries.
            // Import our sample
            importSample(graph, info, entries, sample, chromSizes)
        }
        
        // Now finish up the writing
        graph.finish()
        
        println("VCF imported")
        
    }
    
    /**
     *
     *   Import a sample into the given SequenceGraphBuilder from a VCF, creating
     *   a list of Sides and a list of SequenceGraphEdges, as well as Sites and
     *   Breakpoints for them to be in, and Alleles they contain.
     *
     *   Optionally accepts a chromosome sizes map (of integer lengths by
     *   chromosome name), which, if specified, enables the creation of trailing
     *   telomeres (since we know where to put them)
     *   
     *   Takes the VCF header metadata, an iterator over the VCF's variants, and the
     *   name of the sample to import.
     * 
     */
    def importSample(builder: SequenceGraphBuilder, info: VcfInfo, 
        entries: Iterator[(Variant, List[Metadata.Format], 
        List[List[List[VcfValue]]])], sampleName: String, 
        chromSizes: Option[Map[String, Int]]) {
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
        
        // Where did the last variant end on that contig? It's really 1-bast-
        // the-end, and thus we need to start at 1, since we are using 1-based
        // indexing and the telomere occupies 0.
        var lastEnd: Int = 1
        
        for((variant, formats, samples) <- entries
            if variant.filter == Some(FilterResult.Pass)) {
            // For each variant in the file that passed the filters
            
            // TODO: This whole loop is bad code, working around the lack of a
            // scala "continue" or other easy way to see that, after we've
            // gotten some of the way towards processing this variant, we don't
            // want it after all (and instead we want to complain to the user
            // about its having been there). The Right Way to implement all this
            // is probably as a big for comprehension with lots of guards and
            // code in the for part.
            
            // What contig is it on?
            val contig: String = variant.chromosome match {
                case Left(vcfid) => vcfid.toString
                case Right(string) if string startsWith "chr" => string
                case Right(string) => "chr" + string
            }
            
            // Where does it start in the reference?
            val referenceStart = variant.position
            
            // Where does it end in the reference? This is really 1-bast-the-
            // end, so it's also the next thing's start.
            val referenceEnd = referenceStart + variant.reference.size
            
            // What alleles are available? The reference and all the alternates.
            val alleles = Right(variant.reference) :: variant.alternates
            
            if(contig != lastContig && lastContig != null) {
                // We're ending an old contig and starting a new one.
                
                chromSizes.map { (sizes) =>
                    // We know the chromosome length, so we can add some
                    // trailing telomeres.
                    
                    // How far out do they have to be?
                    val telomereDistance = sizes(lastContig) - lastEnd
                    
                    // Add trailing unphased anchor
                    builder.addAnchor(lastContig, List(0, 1), telomereDistance)
                    
                    // Add trailing telomeres
                    builder.close(lastContig, 0)
                    builder.close(lastContig, 1)
                }
                
                // Create new leading telomeres
                builder.getLastSide(contig, 0)
                builder.getLastSide(contig, 1)
                
                // Record that the leading telomeres are phased, and on this
                // contig, and at the very start.
                lastCallPhased = true
                lastContig = contig
                // The start is at 1 due to 1-based indexing (0, before the
                // beginning, is the telomere, and this is 1-past-the-end)
                lastEnd = 1
            }
            
            // How far are we from the last variant? If this is 0, no Anchors
            // will actually get added.
            val referenceDistance = referenceStart - lastEnd
            
            // This holds the geontype values read for the sample
            val sampleValues = samples(sampleIndex)
            
            // Find the index of the genotype field
            val genotypeIndex = formats indexWhere { (format) =>
                format.id.id == "GT"
            } match {
                // Make it an Option instead of -1 = fail
                case -1 => None
                case somethingElse => Some(somethingElse)
            }
            
            
            // Pull out the genotype string as a String.
            val genotypeString = genotypeIndex.map {(index) =>
                sampleValues(index).head match {
                    case VcfString(value) => value
                    case somethingElse => throw new Exception(
                        "Got a non-string GT value: %s".format(somethingElse))
                }
            }
            
            // Is the string a phased genotype? It is if it hasn't got a "/"
            // in it.
            val genotypePhased = genotypeString.map(_ contains "|")
            
            // Pull out and intify the two allele indices
            val alleleIndices : Option[Seq[Int]] = for(
                genotype <- genotypeString; 
                phased <- genotypePhased;
                result <- {
                    val indices = genotype.split("[|/]").map(_.toInt)
                    if(indices.size != 2) {
                        // Complain about non-diploid places
                        // TODO: Log this.
                        // Skip over them. TODO: work out what they really mean.
                        None
                    } else if(indices.forall(_ == 0) && 
                        phased == lastCallPhased &&
                        !chromSizes.isEmpty) {
                        
                        // They are homozygous reference at this position, and
                        // we don't need this position to change whether we're
                        // phased or not, and we're guaranteed to get an anchor
                        // at the end of the chromosome. So, we can just skip
                        // this position and make it part of the anchor(s) we're
                        // going to add.
                        None
                    } else {
                        // We need to keep going with this position.
                        Some(indices)
                    }
                }
            ) yield result
            
            // At this point, if we don't have None, we know we want this
            // position.
            
            for(phased <- genotypePhased; indices <- alleleIndices) yield {
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
                
                // Report progress
                if(referenceStart % 100 == 0) {
                    // Print progress
                    chromSizes match {
                        case Some(map) =>
                            println("%s:%d: %s (%2.2f%%)".format(contig, 
                                referenceStart, genotype, 
                                (referenceStart:Double)/map(contig) * 100))
                        case None =>
                            println("%s:%d: %s".format(contig, referenceStart, 
                                genotype))
                    }
                }
                
                // Zip allele indices with their phases (confusingly called
                // "indices" also) using zipWithIndex. So we get (allele, phase)
                // tuples.
                val enumeratedIndices = indices.zipWithIndex
                
                
                // Collect by allele index, and pull out phase numbers (so we
                // have all the phases each allele index needs to be stuck in,
                // in a map by allele index).
                val phasesByAlleleIndex = enumeratedIndices.groupBy(_._1)
                    .mapValues(_.map(_._2))
                
                // For phased sites, we look at the phase numbers. For unphased
                // sites, we just look at the lengths.
                
                // For each allele index, get that allele string and add it
                // either to the list of phases, or to all phases with a ploidy
                // given by the list length, depending on phasing status.
                
                // Edge case: for a genotype of just "0" and not "0/0", we'll
                // add an allele just to phase 0.
                phasesByAlleleIndex.foreach { case(alleleIndex, phases) =>
                    // Look up the allele
                    val allele: Either[Either[Breakend, VcfId], String] = 
                        alleles(alleleIndex)
                    
                    // For now we handle only string alleles.
                    val alleleString = allele match {
                        case Right(string) => string
                        case Left(_) => throw new Exception(
                            "Breakends and ID'd alleles not supported")
                    }
                    
                    if(phased) {
                        // Add a different AlleleGroup to each phase
                        phases.map { (phase) =>
                            builder.addAllele(contig, List(phase), 
                                new Allele(alleleString, variant.reference), 
                                referenceEnd - referenceStart)
                        }
                    } else {
                        // Add one AlleleGroup to all of the phases listed, with
                        // a ploidy equal to the list length. We handle the fact
                        // that, at this point, we can't tell the difference
                        // between any phases by having some intervening single
                        // Anchor in both phases between the AlleleGroups for
                        // this variant and anything else.
                        builder.addAllele(contig, phases, 
                            new Allele(alleleString, variant.reference), 
                            referenceEnd - referenceStart)
                    }
                }
                
                // Save information about this site
                lastCallPhased = phased
                lastContig = contig
                lastEnd = referenceEnd
            }
            
        }
        
        // Now we've gone through all variants in the file.
        
        if(lastContig != null) {
            // We had a last contig, so maybe we want trailing telomeres on it?
            chromSizes.map { (sizes) =>
                // We know the chromosome length, so we can add some
                // trailing telomeres.
                
                // How far out do they have to be? The chromosome length is
                // number of actual bases, and we're using 1-based indexing, so
                // we need to add 1 here to get long enough.
                val telomereBases = sizes(lastContig) - lastEnd + 1
                
                // Add trailing unphased anchor
                builder.addAnchor(lastContig, List(0, 1), telomereBases)
                
                // Add trailing telomeres
                builder.close(lastContig, 0)
                builder.close(lastContig, 1)
            }
        }
    }
}





















