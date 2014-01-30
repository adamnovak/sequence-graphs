package edu.ucsc.genome
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

// We want to write files
import java.io._

// We need to work with Avro things
import org.apache.avro.generic.IndexedRecord

// import spark
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD

// import parquet
import parquet.hadoop.{ParquetOutputFormat, ParquetInputFormat}
import parquet.avro.{AvroParquetOutputFormat, AvroWriteSupport, 
                     AvroReadSupport, AvroParquetWriter}

// We need to make Paths for Parquet output.
import org.apache.hadoop.fs.Path

// And we use a hack to get at the (static) schemas of generic things
import org.apache.avro.Schema

// import hadoop job
import org.apache.hadoop.mapreduce.Job


/**
 * SequenceGraphWriter: something that provides a place to put sequence graph
 * components. Write methods must not be called multiple times on the same
 * object.
 */
trait SequenceGraphWriter {
    /**
     * Write the given AlleleGroup.
     */
    def writeAlleleGroup(alleleGroup: AlleleGroup) : Unit
    
    /**
     * Write the given Adjacency.
     */
    def writeAdjacency(adjacency: Adjacency) : Unit
    
    /**
     * Write the given Anchor.
     */
    def writeAnchor(anchor: Anchor) : Unit
    
    /**
     * Add the given Side to the graph. Called exactly once per Side.
     */
    def writeSide(side: Side) : Unit
    
    /**
     * Finish writing the sequence graph.
     */
    def close() : Unit
}

/**
 * GraphvizSequenceGraphWriter: save a sequence graph on disk, writing to a
 * GraphViz dot file in a streaming fashion.
 *
 * Write it to the given file.
 *
 */
class GraphvizSequenceGraphWriter(file: String) extends SequenceGraphWriter {
    
    // We need to write escaped strings. See
    // <http://stackoverflow.com/a/9914380/402891> and
    // <http://commons.apache.org/proper/commons-
    // lang/apidocs/org/apache/commons/lang3/StringEscapeUtils.html>
    import org.apache.commons.lang.StringEscapeUtils.escapeJava
    
    // Open the file for writing. See
    // <http://www.tutorialspoint.com/scala/scala_file_io.htm>
    val graphWriter = new java.io.PrintWriter(new File(file))
    
    // Write a header
    graphWriter.write("digraph sample {\n")
    
    // Set up edge styles
    graphWriter.write("edge [arrowsize=\"0\"];\n")
    
    // And node styles
    graphWriter.write("node [shape=\"point\"];\n")
    
    // Implementations of what we need to be an EasySequenceGraphBuilder
    
    def writeAlleleGroup(alleleGroup: AlleleGroup) : Unit = {
        // What should we label it?
        val label = (alleleGroup.allele match {
            // Put bases if we have actual bases
            case allele: Allele => allele.bases
            // Put something else if we're referencing an allele by index
            case index: java.lang.Integer => "#%d".format(index)
            // It's something else
            case _ => "<unknown>"
        }) + "\nx%d".format(alleleGroup.ploidy.lower)
        
        // Make an edge for every AlleleGroup.
        // TODO: escape labels
        graphWriter.write(
            "%d -> %d [label=\"%s\",arrowsize=\"1\",color=\"#ff0000\"];\n"
            .format(alleleGroup.edge.left, alleleGroup.edge.right, 
            escapeJava(label)))
    }
    
    def writeAdjacency(adjacency: Adjacency) : Unit = {
        // Make an edge for every Adjacency
        graphWriter.write("%d -> %d;\n".format(
            adjacency.edge.left, adjacency.edge.right))
    }
    
    def writeAnchor(anchor: Anchor) : Unit ={
        // What should we write on the anchor?
        val label = "Anchor x%d".format(anchor.ploidy.lower)
        
        // Make an edge for every Anchor
        graphWriter.write("%d -> %d [label=\"%s\",color=\"#0000ff\"];\n"
            .format(anchor.edge.left, anchor.edge.right, escapeJava(label)))
    }
    
    def writeSide(side: Side) : Unit = {
        // What should we write on the Side?
        val label = "%s:%d-%s".format(side.position.contig, side.position.base, 
            side.position.face match {
                case Face.LEFT => "L"
                case Face.RIGHT => "R"
            })
            
        // Getting the label on the point is tricky. See 
        // <http://marc.info/?l=graphviz-interest&m=126029250410816>
            
        // Write out a node to label the side (because point nodes can't have
        // real labels)
        graphWriter.write("L%d [shape=\"plaintext\",label=\"%s\"];\n"
            .format(side.id, escapeJava(label)))

        // Write out the Side without a label
        graphWriter.write("%d;\n"
            .format(side.id))
            
        // Connect them with an invisible edge
        graphWriter.write(
            "{rank=same %d -> L%d [dir=\"none\",style=\"dotted\"];}\n"
            .format(side.id, side.id))
    }
    
    def close() {
        // Close the graph block
        graphWriter.write("}\n")
        
        // Close the file
        graphWriter.close()
        
        println("Wrote %s".format(file))
    }
}

/**
 * ParquetSequenceGraphWriter: a class that writes a sequence graph to Parquet
 * files in the specified directory on the local filesystem. Does not use Spark,
 * just Parquet Avro directly, doing all writing single-threaded to single
 * files.
 *
 * If we used Spark, we would have to load all our data into one RDD before
 * writing it. We could use Spark Streaming to create RDDs batching sequence
 * graph components as they are generated, and write those, but that solution is
 * a bit heavy-weight when we're doing all the graph building in a single thread
 * anyway.
 *
 * Using Spark would enable us to get free access to the appropriate HDFS
 * configuration.
 *
 * If the given directory exists, it will be deleted and replaced.
 *
 */
class ParquetSequenceGraphWriter(directory: String) 
    extends SequenceGraphWriter {
    
    // Remove the output directory, if it exists.
    val directoryFile = new File(directory)
    if(directoryFile.exists()) {
        // We need to delete the directory and everything inside it.
        // Unfortunately this is hard, so we need to use an external library.
        // See <http://alvinalexander.com/blog/post/java/java-io-faq-how-delete-
        // directory-tree>
        import org.apache.commons.io.FileUtils
        FileUtils.deleteDirectory(directoryFile)
    }
    
    
    // Start up writers to the appropriate places. See
    // <https://github.com/bigdatagenomics/adam/blob/master/adam-
    // cli/src/main/scala/edu/berkeley/cs/amplab/adam/cli/Bam2Adam.scala>
    
    // A writer for Sides, which go in /Sides
    val sideWriter = new AvroParquetWriter[Side](
        new Path(directory + "/Sides/sides.parquet"), Side.SCHEMA$)
        
    // A writer for Adjacencies, which go in /Adjacencies
    val adjacencyWriter = new AvroParquetWriter[Adjacency](
        new Path(directory + "/Adjacencies/adjacencies.parquet"),
        Adjacency.SCHEMA$)
        
    // A writer for AlleleGroups, which go in /AlleleGroups
    val alleleGroupWriter = new AvroParquetWriter[AlleleGroup](
        new Path(directory + "/AlleleGroups/AlleleGroups.parquet"),
        AlleleGroup.SCHEMA$)
        
    // A writer for Anchors, which go in /Anchors
    val anchorWriter = new AvroParquetWriter[Anchor](
        new Path(directory + "/Anchors/anchors.parquet"),
        Anchor.SCHEMA$)
    
    def writeAlleleGroup(alleleGroup: AlleleGroup) : Unit = {
        alleleGroupWriter.write(alleleGroup)
    }
    
    def writeAdjacency(adjacency: Adjacency) : Unit = {
        adjacencyWriter.write(adjacency)
    }
    
    def writeAnchor(anchor: Anchor) : Unit = {
        anchorWriter.write(anchor)
    }
    
    def writeSide(side: Side) : Unit = {
        sideWriter.write(side)
    }
    
    def close() {
        // Properly close all the parquet writers, so data gets written.
        sideWriter.close()
        adjacencyWriter.close()
        alleleGroupWriter.close()
        anchorWriter.close()
    }
}

/**
 * SparkParquetSequenceGraphWriter: a class that writes a sequence graph into
 * Parquet files in the specified directory anywhere that Spark can access,
 * using Spark.
 *
 * Requires a SparkContext.
 *
 * Optionally takes a limit (number of objects of each type to accumulate on the
 * master before making a new RDD), and a number of partitions for the RDDs to
 * use.
 *
 * Keeps a limited number of items in memory on the driver. Adds everything in
 * to growing Spark RDDs when the master generates too many things. At the end,
 * saves the Spark RDDs.
 *
 */
class SparkParquetSequenceGraphWriter(directory: String, sc: SparkContext, 
    limit: Int = 10000, partitions: Int = 1) extends SequenceGraphWriter {
    
    // Keep Lists of Sides, Adjacencies, AlleleGroups, and Anchors, when we
    // don't yet have enough to justify creating new RDDs.
    var sides: List[Side] = Nil
    var adjacencies: List[Adjacency] = Nil
    var alleleGroups: List[AlleleGroup] = Nil
    var anchors: List[Anchor] = Nil
    
    // Keep RDDs of each
    var sidesRDD: Option[RDD[Side]] = None
    var adjacenciesRDD: Option[RDD[Adjacency]] = None
    var alleleGroupsRDD: Option[RDD[AlleleGroup]] = None
    var anchorsRDD: Option[RDD[Anchor]] = None
    
    /**
     * Check each type of object to see if we have limit of them or more, and,
     * if so, move the objects from the list to the RDD.
     */
    protected def updateRDDs(itemLimit: Int) = {
    
        if(sides.size >= itemLimit) {
            // Too many Sides. Update the Sides RDD
            
            // Make a new RDD of just the Sides that need to be added.
            val newSidesRDD = sc.makeRDD(sides, partitions)
            
            // Union it in with the old one
            sidesRDD = Some(sidesRDD match {
                case Some(rdd) => rdd.union(newSidesRDD)
                case None => newSidesRDD
            })
            
            // Clear the sides list
            sides = Nil
            
        }
        
        if(adjacencies.size >= itemLimit) {
            // Too many Adjacencies. Update the Sides RDD
            
            // Make a new RDD of just the Adjacencies that need to be added.
            val newAdjacenciesRDD = sc.makeRDD(adjacencies, partitions)
            
            // Union it in with the old one
            adjacenciesRDD = Some(adjacenciesRDD match {
                case Some(rdd) => rdd.union(newAdjacenciesRDD)
                case None => newAdjacenciesRDD
            })
            
            // Clear the adjacencies list
            adjacencies = Nil
            
        }
        
        if(alleleGroups.size >= itemLimit) {
            // Too many AlleleGroups. Update the AlleleGroups RDD
            
            // Make a new RDD of just the AlleleGroups that need to be added.
            val newAlleleGroupsRDD = sc.makeRDD(alleleGroups, partitions)
            
            // Union it in with the old one
            alleleGroupsRDD = Some(alleleGroupsRDD match {
                case Some(rdd) => rdd.union(newAlleleGroupsRDD)
                case None => newAlleleGroupsRDD
            })
            
            // Clear the alleleGroups list
            alleleGroups = Nil
            
        }
        
        if(anchors.size >= itemLimit) {
            // Too many Anchors. Update the Anchors RDD
            
            // Make a new RDD of just the Anchors that need to be added.
            val newAnchorsRDD = sc.makeRDD(anchors, partitions)
            
            // Union it in with the old one
            anchorsRDD = Some(anchorsRDD match {
                case Some(rdd) => rdd.union(newAnchorsRDD)
                case None => newAnchorsRDD
            })
            
            // Clear the anchors list
            anchors = Nil
            
        }
    
    }
    
    /**
     * Given an RDD of Avro records of type RecordType, writes them to a
     * Parquet file in the specified directory. Any other Parquet data in that
     * directory will be overwritten.
     *
     * The caller must have set up Kryo serialization with
     * `SequenceGraphKryoProperties.setupContextProperties()`, so that Avro
     * records can be efficiently serialized to send them to the Spark workers.
     *
     */
    def writeRDDToParquet[RecordType <: IndexedRecord](things: RDD[RecordType],
        directory: String)(implicit m: Manifest[RecordType]) = {
        
        // Every time we switch Avro schemas, we need a new Job. So we create
        // our own.
        val job = new Job()
                
        // Set up Parquet to write using Avro for this job
        ParquetOutputFormat.setWriteSupportClass(job, classOf[AvroWriteSupport])
        
        // Go get the static SCHEMA$ for the Avro record type. We can't do it
        // the normal way (RecordType.SCHEMA$) because Scala keeps static things
        // in companion objects with the same names as the types they really
        // belong to, and type parameters don't automatically bring them along.
        
        // First we need the class of the record, which we need anyway for doing
        // the writing.
        val recordClass = m.erasure
        
        println("Writing out an RDD of %s".format(
            recordClass.getCanonicalName))
        
        // Get the schema Field object    
        val recordSchemaField = recordClass.getDeclaredField("SCHEMA$")
        // Get the static value of this field. Argument is ignored for static
        // fields, so pass null. See <http://docs.oracle.com/javase/7/docs/api/j
        // ava/lang/reflect/Field.html#get%28java.lang.Object%29>
        val recordSchema = recordSchemaField.get(null).asInstanceOf[Schema]
        
        // Set the Avro schema to use for this job
        AvroParquetOutputFormat.setSchema(job, recordSchema)
        
        // Make a PairRDD of serializeable-wrapped Avro records, with null keys.
        val pairRDD = things.map(thing => (null, thing))
        
        // Save the PairRDD to a Parquet file in our output directory. The keys
        // are void, the values are RecordTypes, the output format is a
        // ParquetOutputFormat for RecordTypes, and the configuration of the job
        // we set up is used.
        pairRDD.saveAsNewAPIHadoopFile(directory, classOf[Void], 
            recordClass, classOf[ParquetOutputFormat[RecordType]],
            job.getConfiguration)
                    
    }
    
    // Implementations of what we need to be a SequenceGraphWriter. Just load
    // everything into the appropriate list, and check to see if we need to
    // update the RDDs.
    
    def writeAlleleGroup(alleleGroup: AlleleGroup) : Unit = {
        alleleGroups = alleleGroup :: alleleGroups
        updateRDDs(limit)
    }
    
    def writeAdjacency(adjacency: Adjacency) : Unit = {
        adjacencies = adjacency :: adjacencies
        updateRDDs(limit)
    }
    
    def writeAnchor(anchor: Anchor) : Unit = {
        anchors = anchor :: anchors
        updateRDDs(limit)
    }
    
    def writeSide(side: Side) : Unit = {
        sides = side :: sides
        updateRDDs(limit)
    }
    
    def close() {
        println("Making sure everything is in Spark memory:")
    
        // Make sure everything is in the RDDs
        updateRDDs(0)
        
        // Save everything to the appropriate Parquet directories
        writeRDDToParquet(sidesRDD.get, directory + "/Sides")
        writeRDDToParquet(adjacenciesRDD.get, directory + "/Adjacencies")
        writeRDDToParquet(alleleGroupsRDD.get, directory + "/AlleleGroups")
        writeRDDToParquet(anchorsRDD.get, directory + "/Anchors")

    }
}

/**
 *
 * InMemorySequenceGraphWriter: a class to keep a sequence graph for a genome in
 * memory. Can let the parts of the graph be interrogated.
 *
 */
class InMemorySequenceGraphWriter extends SequenceGraphWriter {

    // We want to write GraphViz graphs through an API.
    import org.kohsuke.graphviz._

    // This holds all the Sides we have created, by ID
    private val sides = HashMap.empty[Long, Side]
    
    // This holds all the AlleleGroups we have created.
    private val alleleGroups = HashSet.empty[AlleleGroup]
    
    // This holds all the Adjacencies we have created.
    private val adjacencies = HashSet.empty[Adjacency]
    
    // This holds all the Anchors we have created.
    private val anchors = HashSet.empty[Anchor]
    
    def writeAlleleGroup(alleleGroup: AlleleGroup) : Unit = 
        alleleGroups += alleleGroup
    
    def writeAdjacency(adjacency: Adjacency) : Unit = 
        adjacencies += adjacency
    
    def writeAnchor(anchor: Anchor) : Unit =
        anchors += anchor
    
    def writeSide(side: Side) : Unit = {
        // Put the Side in the HashMap
        sides(side.id) = side
    }
    
    def close() = {}
    
    // Extra methods we support
    
    /**
     * Get the AlleleGroup attached to the Side with the given ID, or None if no
     * such AlleleGroup exists.
     */
    def getAlleleGroup(sideId: Long) : Option[AlleleGroup] = {
        
        // Find an AlleleGroup that has this ID as one of its edge's ends.
        // There ought to be at most one.
        alleleGroups.find({ (alleleGroup) =>
            (alleleGroup.edge.left == sideId) || 
            (alleleGroup.edge.right == sideId)
        })
    }
    
    /**
     * Write out the graph we have built to the given GraphViz dot file.
     */
    def writeDotFile(file: String) {
        // Make a new graph
        val graph = new org.kohsuke.graphviz.Graph()
        
        // Set up edge styles
        graph.edgeWith(new Style().attr("arrowsize", "0"))
        
        // And node styles
        graph.nodeWith(new Style().attr("shape", "point").attr("label", ""))
        
        // Make a Node for every Side, organized by ID, and add to the graph
        val nodes = sides.mapValues({ (side) =>
            val node = new Node().id(side.id.toString)
            graph.node(node)
            node
        })
        
        alleleGroups map { (alleleGroup : AlleleGroup) =>
            
            // What should we label it?
            val label = (alleleGroup.allele match {
                // Put bases if we have actual bases
                case allele: Allele => allele.bases
                // Put something else if we're referencing an allele by index
                case index: java.lang.Integer => "#%d".format(index)
                // It's something else
                case _ => "<unknown>"
            }) + "\nx%d".format(alleleGroup.ploidy.lower)
            
            // Make an edge for every AlleleGroup.
            graph.edge(nodes(alleleGroup.edge.left),
                nodes(alleleGroup.edge.right), new Style()
                    // Label it with the bases
                    .attr("label", label)
                    // Give it an arrow head
                    .attr("arrowsize", "1")
                    .attr("color", "#ff0000"))
        }
        
        adjacencies map { (adjacency) =>
            // Make an edge for every Adjacency
            graph.edge(nodes(adjacency.edge.left),
                nodes(adjacency.edge.right))
        }
        
        anchors map { (anchor) =>
            // What are the two Sides of the anchor?
            val leftPos = sides(anchor.edge.left).position
            val rightPos = sides(anchor.edge.right).position
            
            // What should we write on the anchor?
            val label = if(leftPos.contig == rightPos.contig) {
                // The anchor is properly on a single contig, so we can
                // determine its length.
                val length = rightPos.base - leftPos.base
                "%sbp anchor x%d".format(length, anchor.ploidy.lower)
            } else {
                // We have no idea how long the anchor is
                "anchor across contigs x%d".format(anchor.ploidy.lower)
            }
            
            // Make an edge for every Anchor
            graph.edge(nodes(anchor.edge.left),
                nodes(anchor.edge.right), new Style()
                    .attr("label", label)
                    .attr("color", "#0000ff"))
        }
        
        // Write the graph to the file
        graph.writeTo(new FileOutputStream(file))
        
        println("Wrote %s".format(file))
        
    }
    
    /**
     * Given an iterable of Avro records of type RecordType, writes them to a
     * Parquet file in the specified directory. Any other Parquet data in that
     * directory will be overwritten.
     *
     * The caller must have set up Kryo serialization with
     * `SequenceGraphKryoProperties.setupContextProperties()`, so that Avro
     * records can be efficiently serialized to send them to the Spark workers.
     *
     */
    def writeCollectionToParquet[RecordType <: IndexedRecord](
        things: Iterable[RecordType], directory: String, sc: SparkContext)
        (implicit m: Manifest[RecordType]) = {
        
        // Every timne we switch Avro schemas, we need a new Job. So we create
        // our own.
        val job = new Job()
        
        try {
            // That SparkContext needs to have the slf4j StaticLoggerBinder we
            // are using on its classpath, or it will complain to standard error
            // about not having that. Since it doesn't inherit our full
            // classpath by default, we need to add the JAR we would load
            // org.slf4j.impl.StaticLoggerBinder from. This is complicated by
            // the fact that that class my not necessarily have actually been
            // provided by the application developer, so we can't just load it
            // and look where it came from.
            
            // Get the class object for the binder we want to send out, via
            // reflection. If it does not exist, an exception will be produced.
            val binderClass = Class.forName("org.slf4j.impl.StaticLoggerBinder")
            
            // Get the .jar it lives in. If it somehow doesn't live in a .jar,
            // this will error out. See
            // <http://stackoverflow.com/q/1983839/402891>
            val binderJar = binderClass.getProtectionDomain.getCodeSource
                .getLocation
            
            // Add the .jar to the Spark context
            sc.addJar(binderJar.toString)
            
        } catch {
            case e: Exception => {
                println("""Failed to find an slf4j StaticLoggerBinder .jar to 
                ship out with Spark jobs. Are you sure you have set up slf4j 
                correctly in your application?""".stripMargin)
                println(e.getMessage)
            }
        }
        
        // TODO: The above shipping out of the log binder .jar *should* work,
        // but does not appear to: the angry message still gets printed by the
        // worker. Figure out why this is and fix it.
                
        // Set up Parquet to write using Avro for this job
        ParquetOutputFormat.setWriteSupportClass(job, classOf[AvroWriteSupport])
        
        // Go get the static SCHEMA$ for the Avro record type. We can't do it
        // the normal way (RecordType.SCHEMA$) because Scala keeps static things
        // in companion objects with the same names as the types they really
        // belong to, and type parameters don't automatically bring them along.
        
        // First we need the class of the record, which we need anyway for doing
        // the writing.
        val recordClass = m.erasure
        
        println("Writing out a collection of %s".format(
            recordClass.getCanonicalName))
        
        // Get the schema Field object    
        val recordSchemaField = recordClass.getDeclaredField("SCHEMA$")
        // Get the static value of this field. Argument is ignored for static
        // fields, so pass null. See <http://docs.oracle.com/javase/7/docs/api/j
        // ava/lang/reflect/Field.html#get%28java.lang.Object%29>
        val recordSchema = recordSchemaField.get(null).asInstanceOf[Schema]
        
        // Set the Avro schema to use for this job
        AvroParquetOutputFormat.setSchema(job, recordSchema)
        
        // Make an RDD of serializeable-wrapped Avro records, with null keys.
        // Put it all into one partition. We have to make sure we make a seq of
        // the values and map on that; if we are not careful with the types, the
        // keys from the HashMap can end up being our RDD keys, and we for some
        // reason need to definitely have a PairRDD (i.e. an RDD of key, value
        // pairs) that definitely has null keys.
        val rdd = sc.makeRDD(things.toSeq map { r => 
            (null, r)
        }, 1)
        
        // Save the RDD to a Parquet file in our output directory. The keys are
        // void, the values are RecordTypes, the output format is a
        // ParquetOutputFormat for RecordTypes, and the configuration of the job
        // we set up is used.
        rdd.saveAsNewAPIHadoopFile(directory, classOf[Void], 
            recordClass, classOf[ParquetOutputFormat[RecordType]],
            job.getConfiguration)
                    
    }
    
    /**
     * Write Parquet Avro files with all of the parts of the graph to the
     * specified directory.
     */
    def writeParquetFiles(directory: String, sc: SparkContext) = {
        // Put the sides under "/Sides"
        writeCollectionToParquet(sides.values, directory + "/Sides", sc)
        // AlleleGroups go under /AlleleGroups
        writeCollectionToParquet(alleleGroups, directory + "/AlleleGroups", sc)
        // Adjacencies go under /Adjacencies
        writeCollectionToParquet(adjacencies, directory + "/Adjacencies", sc)
        // Anchors go under /Anchors
        writeCollectionToParquet(anchors, directory + "/Anchors", sc)
    }
    
    
}

