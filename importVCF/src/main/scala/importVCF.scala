package edu.ucsc.genome.ImportVCF

import scala.collection.mutable.HashMap

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import edu.ucsc.genome._

import java.io.File
import java.nio.file.Paths

// We want to read VCF files
import ca.innovativemedicine.vcf._
import ca.innovativemedicine.vcf.parsers._

// We want to read VCF files in parallel with Hadoop things
import fi.tkk.ics.hadoop.bam.VCFInputFormat
import fi.tkk.ics.hadoop.bam.VariantContextWritable
import org.apache.hadoop.io.LongWritable
import org.broadinstitute.variant.variantcontext.VariantContext


// We want to be able to loop over Java iterators: load a bunch of conversions.
// See <http://stackoverflow.com/a/1625252/402891>
import scala.collection.JavaConversions._

import scala.math._
import scala.util.Sorting

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

object ImportVCF {
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
            val vcfFile = trailArg[String](required = true,
                descr = "VCF file to open")
            // What sample should we import?
            val sampleName = trailArg[String](required = true,
                descr = "Sample to import")
                
            val cluster = opt[String](default = Some("local"),
                descr = "Run against this cluster URL") 
            
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
        
        // Set up logging
        
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
        
        // Set up Spark.
        
        // Set up serialization stuff for Spark so it can efficiently exchange
        // our Avro records.
        SequenceGraphKryoProperties.setupContextProperties()

        // The first thing we need is a Spark context. We would like to be able
        // to make one against any Spark URL: either "local" or soemthing like
        // "mesos://wherever.biz:1234". We need to feed it all the jars it needs
        // to run our code. Fortunately, they all live next to our jar, if the
        // sbt native packager has worked correctly.
        
        // What File is the jar that this class is from?
        val jarFile = new File(getClass.getProtectionDomain
            .getCodeSource.getLocation.toURI)
        
        // What files are in that directory (should all be .jars)? Make a list
        // of their string paths.
        val jarsToSend = jarFile.getParentFile.listFiles.map(_.toString).toSeq
            
        // Set up Spark, giving it the appropriate cluster URL, the SPARK_HOME
        // environment variable (which must be set!), and the list of jars we
        // have worked out.
        println("Initializing Spark")
        val sc = new SparkContext(opts.cluster.get.get, "importVCF", 
            System.getenv("SPARK_HOME"), jarsToSend)
        println("Spark initialized")
        val job = new Job(sc.hadoopConfiguration)
        
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
        
        
        
        // Make a new SequenceGraphWriter to save its graph.
        if(opts.parquetDir.get isDefined) {
            // The user wants to write Parquet. TODO: what if they also asked
            // for a dot file?
            println("Writing to Parquet")
            
            // Import the sample in parallel
            importSampleInParallel(opts.vcfFile.get.get, sample, sc)
                // Write the resulting parts to the directory specified
                .writeToParquet(opts.parquetDir.get.get)
            
            new SparkParquetSequenceGraphWriter(opts.parquetDir.get.get, sc)
        } else if (opts.dotFile.get isDefined) {
            // The user wants to write GraphViz
            println("Writing to Graphviz")
            
            // Make the GraphViz writer
            val writer = new GraphvizSequenceGraphWriter(opts.dotFile.get.get)
            
            // Make the graph builder that writes to the writer
            val graph = new EasySequenceGraphBuilder(sample, "reference", 
                writer)
            
            // Get the actual File or die trying, and then make VcfParser parse
            // that. We need to invoke the VcfParser functor first for some
            // reason.
            VcfParser().parseFile(new File(opts.vcfFile.get.get), false) { 
                (info, entries) =>
                // We get the VCF metadata and an iterator of VCF entries.
                // Import our sample
                importSample(graph, info, entries, sample, chromSizes)
            }
            
            // Now finish up the writing
            graph.finish()
            
        } else {
            // The user doesn't know what they're doing
            throw new Exception("Must specify --dot-file or --parquet-dir")
        }
        
        println("VCF imported")
        
    }
    
    /**
     * Given a VCF Format or Info metadata field list, a list of values for
     * those fields, and the name of the desired field, find the appropriate
     * metadata item for the field, and get the corresponding value(s). If no
     * field with the given name is found, returns None.
     */
    def getFieldValues(formats: Seq[Metadata.HasID], 
        sampleValues: Seq[List[VcfValue]], fieldName: String):
        Option[List[VcfValue]] = {
        
        for(
            // Find the index of the field
            index <- formats indexWhere { (format) =>
                format.id.id == fieldName
            } match {
                // Make it an Option instead of -1 = fail
                case -1 => None
                case somethingElse => Some(somethingElse)
            };
            
            // Pull out the value list
            result <- Some(sampleValues(index))
        ) 
        // Give back the value list. Will give back None if we don't find the
        // value list.
        yield result
    }
    
    /**
     * Given a map from VCF Format or Info metadata fields to lists of values,
     * look up the appropriate value list by ID.
     */
    def getFieldValues(map: Map[_ <: Metadata.HasID, List[VcfValue]], 
        fieldName: String): Option[List[VcfValue]] = {
        
        // Pull apart the metadatas and value lists. It doesn't really help to
        // have an efficient lookup by things we can't construct.
        val (metadatas, values) = map.unzip
        
        // Go look up by list traversal.
        getFieldValues(metadatas.toSeq, values.toSeq, fieldName)
        
    }
    
    /**
     * Function to pair up each element with a copy of its successor. RDD items
     * must be sequentially numbered.
     */
    def pairUp[V: ClassManifest]
        (rdd: RDD[(Long, V)]): RDD[((Long, V), (Long, V))] = {
        
        // Give everything the key of its successor.
        val keyedBySuccessors = rdd.map { (kvPair) =>
            (kvPair._1 + 1, kvPair._2)
        }
        
        // Join up with the successors, transform, and return. We ought to drop
        // the first and last things.
        keyedBySuccessors.join(rdd).map { (item) =>
            // Unpack the entries
            val (successorKey: Long, (thing: V, successor: V)) = item
            
            // Return the thing we want: a pair of the thing and its successor,
            // with both IDs.
            ((successorKey - 1, thing), (successorKey, successor))
        }
        
    }
    
    /**
     * Import a sample from a VCF in parallel, produsing a SequenceGraph, which
     * wraps Spark RDDs.
     */
    def importSampleInParallel(filename: String, sample: String, 
        sc: SparkContext): SequenceGraph = {
        

        // Get an RDD of VCF records, sorted by position. Use VariantContext
        // and not VariantContextWriteable as our internal record type. Make
        // sure it has the correct number of partitions.
        val records: RDD[(Long, VariantContext)] = sc.newAPIHadoopFile(filename,
           classOf[VCFInputFormat], classOf[LongWritable], 
           classOf[VariantContextWritable], sc.hadoopConfiguration)
           .map(p => (p._1.get, p._2.get))
           
        // How many partitions should we use? Apparently we need to match this
        // when we zip two things.
        val numPartitions = records.partitions.size
           
        // TODO: debug why sortByKey() is returning empty RDDs. For now just
        // assume our input VCF is sorted.

        // How many IDs should we allocate per record? Needed for Sides and
        // Edges.
        val idsPerRecord = 100
        
        // Make an RDD with one long ID for each record. IDs are sequential.
        val idBlocks: RDD[Long] = sc.makeRDD(1L until records.count,
            numPartitions)
            
        // Number all the records sequentially. They can be sorted into this
        // order now, and each can guess the IDs of its neighbors.
        val numberedRecords = idBlocks.zip(records.values)
            
        // Get phasing status for each variant
        val phased: RDD[(Long, Boolean)] = records mapValues { (record) =>
            record.getGenotype(sample).isPhased()
        }

        // Produce all the AlleleGroups and their endpoints, numbered by record.
        val alleleGroups: RDD[(Long, SequenceGraphChunk)] = numberedRecords
            .map { (parts) =>
                val (idBlock: Long, record: VariantContext) = parts
                
                // What is the next free ID in our block?
                var nextID = idBlock * idsPerRecord
                
                // What SequenceGraphChunk are we using to collect our parts?
                var chunk = new SequenceGraphChunk()
            
                // Get the sample genotype
                val genotype = record.getGenotype(sample)
                
                // Get the alleles
                val alleles = genotype.getAlleles()
                
                alleles.map { (allele) =>
                    // Get a Side for the left side
                    val side1 = new Side(nextID, new Position(record.getChr, 
                        record.getStart, Face.LEFT), false)
                    nextID += 1
                    // Get a Side for the right side
                    val side2 = new Side(nextID, new Position(record.getChr, 
                        record.getEnd, Face.RIGHT), false)
                    nextID += 1
                    
                    // Make an AlleleGroup between them
                    val group = new AlleleGroup(
                        new Edge(nextID, side1.id, side2.id), 
                        new Allele(allele.getBaseString, 
                            record.getReference.getBaseString), 
                        new PloidyBounds(1, null, null), 
                        sample)
                    nextID += 1

                    // TODO: check for ID overflow
                        
                    // Add everything in to the SequenceGraphChunk we are
                    // building.
                    chunk = chunk + side1 + side2 + group
                }
                
                // Return the assembled chunk of sequence graph, keyed by the
                // record it came from.
                (idBlock, chunk)
            }
           
        // Look at adjacent pairs of SequenceGraphChunks and add the Anchors and
        // Adjacencies to wire them together if appropriate.
        val connectors: RDD[SequenceGraphChunk] = pairUp(alleleGroups).map {
            (item) =>
            // Unpack the item tuple
            val ((firstID, firstChunk), (secondID, secondChunk)) = item
            
            // TODO: check phasing
            // For now just tie together corresponding AlleleGroups.
            val pairsToLink = firstChunk.alleleGroups.zip(
                secondChunk.alleleGroups)
            
            
            // We need to know where to make IDs from. Put them in the first
            // chunk's namespace, after the largest AlleleGroup edge ID we find
            // (since those are made last above).
            var nextID = firstChunk.alleleGroups.map(_.edge.id).max + 1
            
            // We need a SequenceGraphChunk to build in
            var chunk = new SequenceGraphChunk()
            
            pairsToLink.foreach { (pair) =>
                // Unpack the two AlleleGroups to link
                val (firstGroup: AlleleGroup, secondGroup: AlleleGroup) = pair
                
                // Where does our Anchor need to start after?
                val prevPos = firstChunk.getSide(firstGroup.edge.right)
                    .get
                    .position
                
                // Where does our Anchor need to end before?
                val nextPos = secondChunk.getSide(secondGroup.edge.left)
                    .get
                    .position
                
                // Make two Sides for an Anchor
                val startSide = new Side(nextID, new Position(prevPos.contig, 
                    prevPos.base + 1, Face.LEFT), false)
                nextID += 1
                val endSide = new Side(nextID, new Position(nextPos.contig, 
                    nextPos.base - 1, Face.RIGHT), false)
                nextID += 1
                    
                // Make the Anchor
                val anchor = new Anchor(
                    new Edge(nextID, startSide.id, endSide.id), 
                    new PloidyBounds(1, null, null), 
                    sample)
                nextID += 1
                    
                // Make the two connecting adjacencies
                val leftAdjacency = new Adjacency(
                    new Edge(nextID, firstGroup.edge.right, startSide.id), 
                    new PloidyBounds(1, null, null), 
                    sample)
                nextID += 1
                    
                val rightAdjacency = new Adjacency(
                    new Edge(nextID, endSide.id, secondGroup.edge.left), 
                    new PloidyBounds(1, null, null), 
                    sample)
                nextID += 1
                
                // Put everything in the SequenceGraphChunk
                chunk = (chunk + startSide + endSide + anchor + leftAdjacency 
                    + rightAdjacency)
            }
            
            // Return the chunk we built
            chunk
        }
        
        
        
        // Use parallel reduction to get the two minimal-position Sides on each
        // contig. Sides in the tuple are kept in sorted order. Because of the
        // way we build the graph, these are guaranteed to be the dangling
        // endpoints. TODO: Implement this as a GraphX operation to find all
        // degree-1 nodes instead.
        val minimalSidesByContig: Map[String, Seq[Side]] = alleleGroups
            .values
            // Grab the minimal Sides for each chunk
            .map(_.getMinimalSides(2))
            .reduce { (map1, map2) =>
                val pairs = for(
                    // For each contig
                    key <- map1.keySet ++ map2.keySet;
                    // Get the four Sides for that contig
                    values <- Some(map1.getOrElse(key, Nil) ++ 
                        map2.getOrElse(key, Nil));
                    // Pick the best two
                    selected <- Some(Sorting.stableSort(values).slice(0, 2))
                ) yield (key, selected.toSeq)
                
                pairs.toMap
            }
            
        // Do a similar thing to get the maximal-position Sides
        val maximalSidesByContig: Map[String, Seq[Side]] = alleleGroups
            .values
            .map(_.getMaximalSides(2))
            .reduce { (map1, map2) =>
                val pairs = for(
                    key <- map1.keySet ++ map2.keySet;
                    values <- Some(map1.getOrElse(key, Nil) ++ 
                        map2.getOrElse(key, Nil));
                    selected <- Some(Sorting.stableSort(values)
                        // We flip around the sorted sides and pick the
                        // furthest-along two instead.
                        .reverse
                        .slice(0, 2))
                ) yield (key, selected.toSeq)
                
                pairs.toMap
            }
        
        // Now prepend/append leading/trailing Anchors and telomeres.
        
        // Where should we put them?
        var endParts = new SequenceGraphChunk()
        
        // What IDs should we use? Start after the whole range we already used.
        var nextID = records.count * idsPerRecord
        
        // We know both sides maps have the same keys, because if there's a
        // minimal element there must be a maximal one.
        for(
            contig <- minimalSidesByContig.keySet;
            minimalSides <- minimalSidesByContig.get(contig);
            maximalSides <- maximalSidesByContig.get(contig)
        ) {
            // Prepend Anchor and telomere to minimal sides
            for(side <- minimalSides) {
                // Make the Anchor sides.
                // Left starts at base 1
                var anchorLeft = new Side(nextID, 
                    new Position(contig, 1, Face.LEFT), false)
                nextID += 1
                
                // Right ends at the base before the leftmost Side already
                // there.
                var anchorRight = new Side(nextID, 
                    new Position(contig, side.position.base - 1, Face.RIGHT), 
                    false)
                nextID += 1
                
                // Add a Telomere side
                var telomere = new Side(nextID, 
                    new Position(contig, 0, Face.RIGHT), false)
                nextID += 1
                
                // Add the Anchor
                var anchor = new Anchor(
                    new Edge(nextID, anchorLeft.id, anchorRight.id), 
                    new PloidyBounds(1, null, null), 
                    sample)
                nextID += 1
                
                // Add the Adjacency to the rest of the graph
                var adjacencyToGraph = new Adjacency(
                    new Edge(nextID, anchorRight.id, side.id), 
                    new PloidyBounds(1, null, null), 
                    sample)
                nextID += 1
                
                // Add the Adjacency to the telomere
                var adjacencyToTelomere = new Adjacency(
                    new Edge(nextID, telomere.id, anchorLeft.id), 
                    new PloidyBounds(1, null, null), 
                    sample)
                nextID += 1
                
                // Put everything in the chunk
                endParts = (endParts + anchorLeft + anchorRight + telomere + 
                    anchor + adjacencyToGraph + adjacencyToTelomere)
            }
            
            // TODO: Add trailing telomeres for maximalSides according to the
            // chromSizes.
            
        }
        
        
            
        // Aggregate into a SequenceGraph.
        val graph = new SequenceGraph(alleleGroups.values
            .union(connectors)
            .union(sc.parallelize(List(endParts)))
        )
        
        

        // Count up number of Sides made
        println("Sides: %d".format(graph.sides.count))
        
        // Return the finished SequenceGraph
        graph
        
       
    }
    
    /**
     *
     *   Import a sample into the given StructuralSequenceGraphBuilder from a
     *   VCF, creating a list of Sides and a list of SequenceGraphEdges, as well
     *   as Sites and Breakpoints for them to be in, and Alleles they contain.
     *
     *   Optionally accepts a chromosome sizes map (of integer lengths by
     *   chromosome name), which, if specified, enables the creation of trailing
     *   telomeres (since we know where to put them)
     *   
     *   Takes the VCF header metadata, an iterator over the VCF's variants, and
     *   the name of the sample to import.
     * 
     */
    def importSample(builder: StructuralSequenceGraphBuilder, info: VcfInfo, 
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
        
        // Where did the last variant end on that contig? It's really 1-based-
        // the-end, and thus we need to start at 1, since we are using 1-based
        // indexing and the telomere occupies 0.
        var lastEnd: Int = 1
        
        // How far is there a deletion until in each phase of the current
        // contig? Don't make any graph in any phase if there's something in
        // here greater than its start position. TODO: what if variants overlap
        // the ends of deletions?
        val deletedUntil = HashMap.empty[Int, Int]
        
        for(
            // Take appart the variant and see if we really want to process it.
            (variant, formats, samples) <- entries
                if variant.filter == Some(FilterResult.Pass);
            
            // What contig is it on?
            contig <- variant.chromosome match {
                case Left(vcfid) => 
                    Some(vcfid.toString)
                case Right(string) if string startsWith "chr" => 
                    Some(string)
                case Right(string) =>
                    Some("chr" + string)
            };
            
            // Where does it start in the reference?
            referenceStart <- Some(variant.position);
            
            // Where does it end in the reference? This is really 1-past-the-
            // end, so it's also the next thing's start.
            referenceEnd <- Some(referenceStart + variant.reference.size);
            
            // What alleles are available? The reference and all the alternates.
            alleles <- Some(Right(variant.reference) :: variant.alternates);
            
            // This holds the geontype values read for the sample
            sampleValues <- Some(samples(sampleIndex));
            
            // Find the index of the genotype field
            genotypeIndex <- formats indexWhere { (format) =>
                format.id.id == "GT"
            } match {
                // Make it an Option instead of -1 = fail
                case -1 => None
                case somethingElse => Some(somethingElse)
            };
            
            // Get the genotype(s) list for the sample. It should contain only
            // one genotype. TODO: So something special if "GT" isn't found.
            genotypes <- getFieldValues(formats, sampleValues, "GT");
            
            // Pull out the genotype string as a String.
            genotypeString <- genotypes.head match {
                case VcfString(value) =>
                    // We found the genotype string.
                    Some(value)
                case somethingElse => 
                    // Somehow the genotype isn't a string?
                    throw new Exception("Got a non-string GT value: %s".format(
                        somethingElse))
                    None
            };
            
            // Is the string a phased genotype? It is if it hasn't got a "/"
            // in it.
            genotypePhased <- Some(genotypeString contains "|");
            
            // Pull out and intify the two allele indices
            alleleIndices <- {
                val indices = genotypeString.split("[|/]").map(_.toInt)
                if(indices.size != 2) {
                    // Complain about non-diploid places
                    // TODO: Log this.
                    // Skip over them. TODO: work out what they really mean.
                    None
                } else if(indices.forall(_ == 0) && 
                    genotypePhased == lastCallPhased &&
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
            };
            
            // Zip allele indices with their phases (confusingly called
            // "indices" also) using zipWithIndex. So we get (allele, phase)
            // tuples.
            enumeratedIndices <- Some(alleleIndices.zipWithIndex);
            
            // Collect by allele index, and pull out phase numbers (so we have
            // all the phases each allele index needs to be stuck in, in a map
            // by allele index). For phased sites, we look at the phase numbers.
            // For unphased sites, we just look at the number of phases.
            phasesByAlleleIndex <- Some(enumeratedIndices.groupBy(_._1)
                .mapValues(_.map(_._2).filter { (phase) =>
                    // Only let through phases that aren't deleted at the
                    // variant's start point
                    referenceStart > deletedUntil.getOrElse(phase, -1)
                }).filter(_._2.size > 0))
            
        ) {
            // For each variant in the file that passed the filters
            
            if(contig != lastContig && lastContig != null) {
                // We're ending an old contig and starting a new one.
                
                chromSizes.map { (sizes) =>
                    // We know the chromosome length, so we can add some
                    // trailing telomeres.
                    
                    // How far out do they have to be?
                    val telomereDistance = sizes(lastContig) - lastEnd
                    
                    // Add trailing unphased anchor
                    builder.addAnchor(lastContig, List(0, 1), lastEnd, 
                        telomereDistance)
                    
                    // Add trailing telomeres
                    builder.close(lastContig, 0, sizes(lastContig))
                    builder.close(lastContig, 1, sizes(lastContig))
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
            
            // Add the intermediate Anchors between the last variant site and
            // this one. If there's no actual space between the sites, it will
            // be a backwards anchor (end is at a 1-base-earlier position than
            // start)
            if(genotypePhased && lastCallPhased) {
                // We need phased anchors, since both this call and the
                // previous one are phased. TODO: We assume a diploid genotype.
                builder.addAnchor(contig, List(0), lastEnd, referenceDistance)
                builder.addAnchor(contig, List(1), lastEnd, referenceDistance)
                
            } else {
                // We need an unphased anchor.
                builder.addAnchor(contig, List(0, 1), lastEnd,
                    referenceDistance)
            }
            
            
            
            // Report progress
            if(referenceStart % 100 == 0) {
                // Print progress
                chromSizes match {
                    case Some(map) =>
                        println("%s:%d: %s (%2.2f%%)".format(contig, 
                            referenceStart, genotypeString, 
                            (referenceStart:Double)/map(contig) * 100))
                    case None =>
                        println("%s:%d: %s".format(contig, referenceStart, 
                            genotypeString))
                }
            }
            
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
                allele match {
                    case Right(alleleString) =>
                        // We have a normal string replacement allele.
                        if(genotypePhased) {
                            // Add a different AlleleGroup to each phase
                            phases.map { (phase) =>
                                builder.addAllele(contig, List(phase),
                                    referenceStart, 
                                    new Allele(alleleString, variant.reference), 
                                    referenceEnd - referenceStart)
                            }
                        } else {
                            // Add one AlleleGroup to all of the phases listed,
                            // with a ploidy equal to the list length. We handle
                            // the fact that, at this point, we can't tell the
                            // difference between any phases by having some
                            // intervening single Anchor in both phases between
                            // the AlleleGroups for this variant and anything
                            // else.
                            builder.addAllele(contig, phases,
                                referenceStart, 
                                new Allele(alleleString, variant.reference), 
                                referenceEnd - referenceStart)
                        }
                        
                        // These sorts of alleles let you continue after the
                        // reference allele's bases in the reference, if they
                        // aren't alternatives to/the reference for something
                        // fancier..
                        lastEnd = max(lastEnd, referenceEnd)
                        
                    case Left(Right(alleleId)) =>
                        // We have a named allele. Ignore the ID, but look at
                        // the SVTYPE and SVLEN to figure out what to do.
                        
                        // Get the (only) structural variant type.
                        // TODO: Stick list logic into getFieldValues?
                        val variantType = getFieldValues(variant.info,
                            "SVTYPE") match {
                                case Some(VcfString(head) :: Nil) => head
                                case somethingElse => throw new Exception(
                                    "Bad structural variant type: %s"
                                        .format(somethingElse))
                            }
                            
                        // Get the (only) structural variant length
                        val variantLength = getFieldValues(variant.info,
                            "SVLEN") match {
                                case Some(VcfInteger(head) :: Nil) => head
                                case somethingElse => throw new Exception(
                                    "Bad structural variant length: %s"
                                        .format(somethingElse))
                            }
                            
                        // Get the (only) structural variant end position, which
                        // is optional.
                        val variantEnd = getFieldValues(variant.info,
                            "END") match {
                                case Some(VcfInteger(head) :: Nil) => Some(head)
                                case somethingElse => None
                            }
                            
                        println("Structural variant: %s".format(variantType))
                        
                        variantType match {
                            case "INS" =>
                                // Process an insert. SVLEN tells us how many
                                // new bases the insert contains, in addition to
                                // the reference allele.
                                
                                // How many bases do we need? Reference length
                                // plus SVLEN, and all of them are N since we
                                // have no sequence.
                                val bases = "N" * (referenceEnd -
                                    referenceStart + variantLength)
                                
                                if(genotypePhased) {
                                    // Add a different AlleleGroup to each phase
                                    phases.map { (phase) =>
                                        builder.addAllele(contig, List(phase),
                                            referenceStart, 
                                            new Allele(bases, 
                                                variant.reference), 
                                            referenceEnd - referenceStart)
                                    }
                                } else {
                                    // Add one AlleleGroup to all of the phases
                                    // listed, with a ploidy equal to the list
                                    // length.
                                    builder.addAllele(contig, phases,
                                        referenceStart, 
                                        new Allele(bases, variant.reference), 
                                        referenceEnd - referenceStart)
                                }
                                
                                // Inserts let you continue right after the
                                // reference allele.
                                lastEnd = max(lastEnd, referenceEnd)
                                
                            case "DEL" =>
                                // Process a deletion. SVLEN should be negative
                                // and say how many bases relative to the
                                // reference allele were lost. So a reference
                                // allele of A with SVLEN=-1 means the alt
                                // allele is "".
                                
                                // Get the recorded END, or try to reconstruct
                                // it.
                                val endpoint = variantEnd.getOrElse {
                                    (referenceStart + -variantLength)
                                }
                                
                                println("Ends at %d".format(endpoint))
                                
                                // Phasing doesn't really matter here. Add as
                                // deletion up to the endpoint.
                                builder.addDeletion(contig, phases,
                                    referenceStart,
                                    endpoint - referenceStart)
                                    
                                // Now remember where we are deleted until the
                                // endpoint on each phase. TODO: This will cause
                                // trouble if e.g. a sample has a small deletion
                                // and a larger spanning deletion, and they are
                                // unphased and both end up on phase 0.
                                phases.map { (phase) =>
                                   deletedUntil(phase) = endpoint
                                }
                                
                                // Deletions make you continue after the implied
                                // longer reference allele that accounts for
                                // SVLEN or END.
                                lastEnd = max(lastEnd, endpoint)
                        }
                            
                        
                    case Left(Left(breakend)) => throw new Exception(
                        "Found breakend: %s".format(breakend))
                }
                
                
            }
            
            if(phasesByAlleleIndex.size == 0) {
                // We didn't actually have any genotypes (somehow)
            }
            
            // Save information about this site
            lastCallPhased = genotypePhased
            lastContig = contig
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
                builder.addAnchor(lastContig, List(0, 1), lastEnd,
                    telomereBases)
                
                // Add trailing telomeres
                builder.close(lastContig, 0, sizes(lastContig))
                builder.close(lastContig, 1, sizes(lastContig))
            }
        }
    }
}





















