package edu.ucsc.genome.ImportVCF

import scala.collection.mutable.HashMap

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import org.apache.spark.rdd.RDD

import org.apache.spark.graphx
import org.apache.spark.graphx._

import edu.ucsc.genome._

import java.io.File
import java.nio.file.Paths

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

import org.apache.hadoop.mapreduce.Job

/**
 * A fully serializeable approximation of VariantContext that doesn't have any
 * lingering lazy parsing to do. Just copies over all the data we actually use.
 *
 * Shouldn't be secretly storing the original VariantContext as a field.
 *
 * TODO: Replace with VariantContextWritable
 */
class VariantRecord(@transient original: VariantContext) extends Serializable {
    import org.broadinstitute.variant.variantcontext.{Allele, Genotype}
    
    // Steal all the values from the original
    val chr: String = original.getChr
    val start: Long = original.getStart
    val end: Long = original.getEnd
    val alleles: Seq[Allele] = original.getAlleles
    val reference: Allele = original.getReference
    val genotypes: Map[String, Genotype] = {
        // Go through all the samples and get their genotypes, and make our own
        // map.
        original.getSampleNames.map {
            (name: String) => (name, original.getGenotype(name))
        }.toMap
    }
    val filtered: Boolean = original.isFiltered
    
    
    // Replicate a portion of the interface
    def getChr: String = chr
    def getStart: Long = start
    def getEnd: Long = end
    def getAlleles: Seq[Allele] = alleles
    def getReference: Allele = reference
    def getGenotype(sample: String): Genotype = genotypes(sample)
    def isFiltered: Boolean = filtered
    
}


object ImportVCF {
    def main(args: Array[String]) {
        
        // Option parsing code more or less stolen from various Scallop
        // examples.
        val opts = new ScallopConf(args) {
            guessOptionName = true
            banner("""Usage: importVCF [OPTION] {--dot-file dotFile | 
                 --parquet-dir parquetDir | --serial parquetDir} vcfFile sample
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
                descr = "Save Parquet files in this directory (absolute path)")
            
            val version = opt[Boolean](noshort = true, 
                descr = "Print version")
            val help = opt[Boolean](noshort = true, 
                descr = "Show this message")

        } 
        
        // Set up logging
        ConsoleUtil.quietParquet
        
        // Set up Spark.
        
        // Set up serialization stuff for Spark so it can efficiently exchange
        // our Avro records.
        SequenceGraphKryoProperties.setupContextProperties()
        
        // Set the executors to use a more reasonable amount of memory.
        System.setProperty("spark.executor.memory", "25G")
        
        // Set the parallelism level to have enough reducers to not run out of
        // memory
        System.setProperty("spark.default.parallelism", "12")
        
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
        
        // Import the sample in parallel, but only if we actually try to save
        // it.
        lazy val imported = importSampleInParallel(opts.vcfFile.get.get, sample,
            chromSizes, sc)
        
        // Make a new SequenceGraphWriter to save its graph.
        if(opts.parquetDir.get isDefined) {
            // The user wants to write Parquet. TODO: what if they also asked
            // for a dot file?
            println("Will write to Parquet")
            
            
            // Write the resulting parts to the directory specified
            imported.writeToParquet(opts.parquetDir.get.get)
        } else if (opts.dotFile.get isDefined) {
            // The user wants to write GraphViz
            println("Will write to GraphViz")
            
            // Make the GraphViz writer
            val writer = new GraphvizSequenceGraphWriter(opts.dotFile.get.get)
            
            // Write all the graph elements to it in serial.
            imported.write(writer)
            
            // Now finish up the writing
            writer.close()
        } else {
            // The user doesn't know what they're doing
            throw new Exception("Must specify --dot-file or --parquet-dir")
        }
        
        println("VCF imported")
        
    }
    
    /**
     * Run the given function for each pair of adjacent elements in the given
     * RDD. The first element will be involved in an invocation with None as the
     * first argument, and the last argument will be involved as an invocation
     * with None as the second argument.
     *
     * Implementation based on method provided by Imran Rashid: <http://mail-
     * archives.apache.org/mod_mbox/spark-user/201402.mbox/%3CCAO24D
     * %3DTkwfWbRiBrz2dP84u-cVXG7CT9%2BkoK5dfUMBQC%3D6y3nQ%40mail.gmail.com%3E>
     *
     * Note that the type in the RDD must be properly serializeable. If it does
     * something like lazy parsing (like GATK's VariantContext does), this
     * method will break when it tries to send open file handles or whatever
     * from one place to another.
     */
    def withNeighbors[T: ClassManifest, U: ClassManifest](rdd: RDD[T], 
        function: (Option[T], Option[T]) => U): RDD[U] = {
        
        // We re-use our input rdd, so cache it.
        rdd.cache
        
        // Make a map from partition to the first element of that partition.
        // First make an RDD of partition first element by partition index. Ends
        // up being a 1-element-per-partition RDD.
        val firstElements: RDD[(Int, Option[T])] = rdd.mapPartitionsWithIndex { 
            case(index, iterator) =>
                iterator.hasNext match {
                    // We have a first element. Say so.
                    case true => 
                        val element = iterator.next
                        Some((index, Some(element))).iterator
                    // We have no first element (empty partition). Say so.
                    case false => 
                        Some((index, None)).iterator
                }
            case otherThing => throw new Exception("Wrong thing!")
        }
        
        // Now collect that into an Array of first elements (or None) ordered by
        // partition. Probably easier to do it on the master than via an
        // accumulator that the master would need to read anyway.
        val broadcastArray = rdd.context.broadcast(firstElements
                // Collect to master (depends on T being able to actually
                // serialize/deserialize and not just pretending to)
                .collect
                // Sort by partition
                .sortBy(_._1)
                // Get the first value (or None)
                .map(_._2))

        // Now do a scan in each partition, and produce the results for each
        // partition.
        rdd.mapPartitionsWithIndex { case(index, iterator) =>
            // Put a leading None on the first partition, but nothing everywhere
            // else.
            val leadingIterator = index match {
                case 0 => List(None).iterator
                case _ => Iterator.empty
            }
            
            // Find the next value after our partition. It may not exist, in
            // which case we want None. TODO: Can we have a flat version of
            // this?
            val nextValue: Option[Option[T]] = broadcastArray.value
                // Starting at the thing after us
                .drop(index + 1)
                // Find the first Some and grab it
                .find(!_.isEmpty)
                
            // So we have either Some(Some(value)) or None (for no Somes found)
            
            // Make an iterator for the next value
            val trailingIterator: Iterator[Option[T]] = nextValue match {
                // We found something, so give back an iterator for what we
                // found.
                case Some(Some(found)) => List(Some(found)).iterator
                // There are no non-empty partitions after us, so add a trailing
                // None.
                case None => List(None).iterator
                // Otherwise we broke
                case otherwise => 
                    throw new Exception("Next matching thing was %s".format(
                        otherwise))
            }
                
            // Concatenate the leading None (if needed), all out items as
            // Options, and the trailing None (if needed).
            val optionIterator: Iterator[Option[T]] = leadingIterator ++
                iterator.map(Some(_)) ++ 
                trailingIterator
            
            // Do a window-2 scan over the iterator, sliding over by 1 each time
            // (the default). Skip partial windows. See <http://daily-
            // scala.blogspot.com/2009/11/iteratorsliding.html>
            val slidingIterator = optionIterator.sliding(2, 1)
            val groupedIterator = slidingIterator.withPartial(false)
            
            // Run the function and return an iterator of its results.
            groupedIterator.map { (pair: Seq[Option[T]]) =>
                // We have a length-2 Seq of the two things. Just pull out the
                // two items. No possibility of partial windows, since we banned
                // them.
                function(pair(0), pair(1))
            }
        }
        
        
        
    }
    
    /**
     * Import a sample from a VCF in parallel, producing a SequenceGraph, which
     * wraps Spark RDDs.
     *
     * Takes the VCF file name to import, the name of the sample to look at in
     * that file, an optional map from contig names to sizes, and a Spark
     * context.
     *
     */
    def importSampleInParallel(filename: String, sample: String, 
        chromSizes: Option[Map[String, Int]], 
        sc: SparkContext): SequenceGraph = {
        

        // Get an RDD of VCF records, sorted by position. Use VariantRecord and
        // not VariantContext or VariantContextWriteable as our internal record
        // type, since VariantContext does not want to be serialized and we do
        // need to consider records with their neighbors.
        val records: RDD[(Long, VariantRecord)] = sc.newAPIHadoopFile(filename,
           classOf[VCFInputFormat], classOf[LongWritable], 
           classOf[VariantContextWritable], sc.hadoopConfiguration)
           .map(p => (p._1.get, new VariantRecord(p._2.get)))
           
        println("Prepared VCF")
        
        // How many partitions should we use? Apparently we need to match this
        // when we zip two things.
        val numPartitions = records.partitions.size
           
        // Filter down the records, throwing out any that are "filtered"
        val passingRecords = records.filter { (pair) => !(pair._2.isFiltered) }
        
        // Report how many records pass filters, so we can get a handle on how
        // many records are dropped by the next step.
        println("Counting records passing filters...")
        println("Got %d records".format(passingRecords.count))
        
        // Figure out which records are interesting: actually non-reference in
        // some way, or important in gaining/losing phasing.
        val interestingRecords: RDD[(Long, VariantRecord)] = chromSizes match {
            // We know chromosome sizes, so every record can rely on having an
            // anchor after it.
            case Some(sizes) => withNeighbors(passingRecords, {
                // Go through and look at the last record and the current record
                // for every record.
                (previous: Option[(Long, VariantRecord)],
                    current: Option[(Long, VariantRecord)]) =>

                current match {
                    case Some((id, record)) =>
                        // We have the current record. Grab the genotype.
                        val currentGenotype = record.getGenotype(sample)
                        
                        previous match {
                            case Some((id2, record2)) =>
                                // We have the previous one too. Grab it too.
                                val previousGenotype = record2.getGenotype(
                                    sample)
                                
                                if(currentGenotype.isHomRef && 
                                    currentGenotype.isPhased == 
                                    previousGenotype.isPhased) {
                                    
                                    // It's homozygous reference and doesn't
                                    // change away from our previous phasing.
                                    // Discard it.
                                    None
                                    
                                } else {
                                    // It's either variant or needed to change
                                    // phasing. Keep it
                                    Some((id, record))
                                }
                            case None =>
                                // This is the first record. Keep it.
                                // TODO: See when we can drop the first record.
                                Some((id, record))
                        }
                    // We have no current record. Don't produce anything.
                    case None => None
                }
                
            }).flatMap((x: Option[(Long, VariantRecord)]) => x)
            
            // If chromosome sizes aren't known, we won't generate trailing
            // telomere anchors, so we can't let things get subsumed by the
            // anchors that follow them. Just pass everything through.
            case None => passingRecords
        }
        
        // Go count up the number of records passing filters, which we need to
        // know in order to assign IDs.
        println("Counting exact number of records to convert...")
        val recordCount = interestingRecords.count
        println("Got %d records".format(recordCount))
           
        // TODO: debug why sortByKey() is returning empty RDDs. For now just
        // assume our input VCF is sorted.

        // How many IDs should we allocate per record? Needed for Sides and
        // Edges.
        val idsPerRecord = 100
        
        // Make an RDD with one long ID for each record. IDs are sequential.
        val idBlocks: RDD[Long] = sc.makeRDD(0L until recordCount,
            numPartitions)
            
        // Number all the records sequentially. They can be sorted into this
        // order now, and each can guess the IDs of its neighbors.
        val numberedRecords = idBlocks.zip(interestingRecords.values).cache
            
        // Get phasing status for each variant, keyed by sequential ID.
        val phased: RDD[Boolean] = numberedRecords map { case (_, record) => 
            record.getGenotype(sample).isPhased()
        }
        
        // Produce all the AlleleGroups and their endpoints
        val alleleGroups: RDD[SequenceGraphChunk] = numberedRecords
            .map { (parts) =>
                val (idBlock: Long, record: VariantRecord) = parts
                
                // Make a pen to draw Sequence Graph parts with IDs from the
                // given block.
                val pen = new SequenceGraphPen(sample, idBlock * idsPerRecord,
                    idsPerRecord)
                
                // What SequenceGraphChunk are we using to collect our parts?
                var chunk = new SequenceGraphChunk()
            
                // Get the sample genotype
                val genotype = record.getGenotype(sample)
                
                // Get the alleles
                val alleles = genotype.getAlleles()
                
                // Get the contig we're on
                val contig = record.getChr match {
                    // Things that start with "chr" pass as is
                    case string if string startsWith "chr" => string
                    // Everything else gets "chr" prepended
                    case string => "chr" + string
                    // TODO: Handle "random" or "ecoli_whatever" or what have
                    // you.
                };
                
                alleles.map { (vcfAllele) =>
                    // Convert the VCF allele to one of our Alleles, which needs
                    // both sample and reference strings (for export)
                    val allele = new Allele(vcfAllele.getBaseString, 
                            record.getReference.getBaseString)
                    
                    // Get a Side for the left side
                    val side1 = pen.drawSide(contig, record.getStart, Face.LEFT)
                    chunk += side1
                    
                    // Get a Side for the right side
                    val side2 = pen.drawSide(contig, record.getEnd, Face.RIGHT)
                    chunk += side2
                    
                    // Make an AlleleGroup between them
                    chunk += pen.drawAlleleGroup(side1, side2, allele)
                }
                
                // Return the assembled chunk of sequence graph
                chunk
            }.cache
        
        // Get (alleleGroupChunk, phasedFlag) tuples
        val alleleGroupsWithPhasing = alleleGroups.zip(phased)
           
        // Look at adjacent pairs of SequenceGraphChunks and add the Anchors and
        // Adjacencies to wire them together if appropriate.
        // TODO: Adapt the pairUp thing to be able to reject variants after some
        // processing.
        val connectors: RDD[SequenceGraphChunk] = 
            withNeighbors(alleleGroupsWithPhasing, {
             
            (left: Option[(SequenceGraphChunk, Boolean)], 
            right: Option[(SequenceGraphChunk, Boolean)]) =>
            
            // We need a SequenceGraphChunk to build in
            var chunk = new SequenceGraphChunk()
            
            for(
                // Go grab the chunks and their metadata for both sides.
                (firstChunk, firstPhased) <- left;
                (secondChunk, secondPhased) <- right
            ) {
                
                // We need to know where to make IDs from. Put them in the first
                // chunk's namespace, after the largest AlleleGroup edge ID we
                // find (since those are made last above).
                val nextID = firstChunk.alleleGroups.map(_.edge.id).max + 1

                // Make a pen to draw Sequence Graph parts with IDs from the
                // given block. We need to reconstruct where the block ought to
                // end, because we sort of hack about to grab the next ID to
                // use.
                val pen = new SequenceGraphPen(sample, nextID, 
                    idsPerRecord - nextID % idsPerRecord)
                
                
                
                // Find a Side where the first variant ended. TODO: deal with
                // complex variants.
                val prevSideID = firstChunk.alleleGroups.head.edge.right
                
                // Get its Position.
                val prevPos = firstChunk.getSide(prevSideID).get.position
                
                // Find a Side where the next variant started. TODO: deal with
                // complex variants.
                val nextSideID = secondChunk.alleleGroups.head.edge.right
                
                // Get its Position.
                val nextPos = secondChunk.getSide(nextSideID).get.position
                
                if(prevPos.contig == nextPos.contig) {
                    // Need to be physically linked since they're on the same
                    // contig.
                    if(firstPhased && secondPhased) {
                        // Link corresponding pairs of alleles with Anchors
                        
                        // TODO: De-nest into its own function or something,
                        // this is super-indented.
                        
                        val pairsToLink = firstChunk.alleleGroups.zip(
                            secondChunk.alleleGroups)
                            
                        pairsToLink.foreach { (pair) =>
                            // Unpack the two AlleleGroups to link
                            val (firstGroup: AlleleGroup, 
                                secondGroup: AlleleGroup) = pair
                            
                            // We use prevPos and nextPos as found for the two
                            // Sites above. This will work as long as variants
                            // don't start generating anything more interesting
                            // than AlleleGroups with identical reference
                            // positions.
                            
                            // Make two Sides for an Anchor
                            val startSide = pen.drawSide(prevPos.contig, 
                                prevPos.base + 1, Face.LEFT)
                            chunk += startSide
                            
                            val endSide = pen.drawSide(nextPos.contig, 
                                nextPos.base - 1, Face.RIGHT)
                            chunk += endSide
                                
                            // Make the Anchor
                            chunk += pen.drawAnchor(startSide, endSide)
                                
                            // Make the two connecting adjacencies
                            chunk += pen.drawAdjacency(firstGroup.edge.right,
                                startSide)
                            chunk += pen.drawAdjacency(endSide,
                                secondGroup.edge.left)
                        }
                        
                    } else {
                        // Lose phasing by linking through a single Anchor
                        
                        // Make two Sides for an Anchor
                        val startSide = pen.drawSide(prevPos.contig, 
                            prevPos.base + 1, Face.LEFT)
                        chunk += startSide
                        
                        val endSide = pen.drawSide(nextPos.contig, 
                            nextPos.base - 1, Face.RIGHT)
                        chunk += endSide
                        
                        // Make the Anchor
                        chunk += pen.drawAnchor(startSide, endSide)
                        
                        firstChunk.alleleGroups.foreach { (alleleGroup) =>
                            // Attach the end of each of the AlleleGroups on the
                            // left
                            chunk += pen.drawAdjacency(alleleGroup.edge.right,
                                startSide)
                        }
                        
                        secondChunk.alleleGroups.foreach { (alleleGroup) =>
                            // Attach the start of each of the AlleleGroups on
                            // the left
                            chunk += pen.drawAdjacency(endSide,
                                alleleGroup.edge.left)
                        }
                        
                        // TODO: Unify this with the phased case, somehow. Or
                        // write utility methods (leftmost/rightmost Position
                        // from a chunk?).
                        
                        // Now we have hooked the two chunks together, losing
                        // phasing.
                        
                    }
                }
            
            }
            
            // Return the chunk (which has nothing in it if we didn't feel the
            // need to build anything, because we are, for example, on an end).
            chunk
            
        })
        
        
        // Use parallel reduction to get the two minimal-position Sides on each
        // contig. Sides in the tuple are kept in sorted order. Because of the
        // way we build the graph, these are guaranteed to be the dangling
        // endpoints. TODO: Implement this as a GraphX operation to find all
        // degree-1 nodes instead.
        val minimalSidesByContig: Map[String, Seq[Side]] = alleleGroups
            // Grab the minimal Sides for each chunk
            .map(_.getMinimalSides(2))
            .fold(Map.empty) { (map1, map2) =>
                
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
            .map(_.getMaximalSides(2))
            .fold(Map.empty) { (map1, map2) =>
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
        val nextID = recordCount * idsPerRecord
        
        // Make a pen to draw the stuff on the ends of the chromosomes, in
        // serial. It can go as far as it wants in ID space.
        val pen = new SequenceGraphPen(sample, nextID)
        
        // We know both sides maps have the same keys, because if there's a
        // minimal element there must be a maximal one. There may be neither,
        // but this for comprehension should handle that.
        for(
            contig <- minimalSidesByContig.keySet;
            minimalSides <- minimalSidesByContig.get(contig);
            maximalSides <- maximalSidesByContig.get(contig)
        ) {
            for(side <- minimalSides) {
                // Prepend Anchor and telomere to minimal Sides.
            
                // Make the Anchor sides.
                // Left starts at base 1
                val anchorLeft = pen.drawSide(contig, 1, Face.LEFT)
                endParts += anchorLeft
                
                // Right ends at the base before the leftmost Side already
                // there.
                val anchorRight = pen.drawSide(contig, side.position.base - 1,
                    Face.RIGHT)
                endParts += anchorRight
                
                // Add a Telomere side
                val telomere = pen.drawSide(contig, 0, Face.RIGHT)
                endParts += telomere
                
                // Add the Anchor
                endParts += pen.drawAnchor(anchorLeft, anchorRight)
                
                // Add the Adjacency to the rest of the graph
                endParts += pen.drawAdjacency(anchorRight, side)
                
                // Add the Adjacency to the telomere
                endParts += pen.drawAdjacency(telomere, anchorLeft)
            }
            
            for(sizes <- chromSizes; side <- maximalSides) {
                // Since we have the chromosome sizes, do trailing telomeres on
                // maximal Sides as well.
                
                // Make the Anchor sides.
                // Left starts at base after the rightmost Side already there.
                val anchorLeft = pen.drawSide(contig, side.position.base + 1,
                    Face.LEFT)
                endParts += anchorLeft
                
                // Right ends at the last base
                val anchorRight = pen.drawSide(contig, sizes(contig),
                    Face.RIGHT)
                endParts += anchorRight
                
                // Add a Telomere side
                val telomere = pen.drawSide(contig, sizes(contig) + 1,
                    Face.LEFT)
                endParts += telomere
                
                // Add the Anchor
                endParts += pen.drawAnchor(anchorLeft, anchorRight)
                
                // Add the Adjacency to the rest of the graph
                endParts += pen.drawAdjacency(side, anchorLeft)
                
                // Add the Adjacency to the telomere
                endParts += pen.drawAdjacency(anchorRight, telomere)
                
            }
            
        }
        
        // Aggregate into a SequenceGraph.
        val graph = new SequenceGraph(alleleGroups
            .union(connectors)
            .union(sc.parallelize(List(endParts)))
        )

        println("Constructed final SequenceGraph RDDs")
        
        // Return the finished SequenceGraph
        graph
    }
    
}





















