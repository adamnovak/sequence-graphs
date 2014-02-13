package edu.ucsc.genome

import java.io.{ByteArrayOutputStream, ByteArrayInputStream, DataOutputStream,
    DataInputStream}

import scala.collection.JavaConversions._

import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.graphx._
import org.apache.avro.generic.IndexedRecord

// Import parquet
import parquet.hadoop.{ParquetOutputFormat, ParquetInputFormat}
import parquet.avro.{AvroParquetOutputFormat, AvroWriteSupport, AvroReadSupport,
    AvroParquetWriter}
import parquet.hadoop.util.{ContextUtil, ConfigurationUtil}

// We need to make Paths for Parquet output.
import org.apache.hadoop.fs.Path

// And we use a hack to get at the (static) schemas of generic things
import org.apache.avro.Schema

// And more hacks to get things into the right type when Parquet Avro ignores
// our classpath.
import org.apache.avro.io.{EncoderFactory, DecoderFactory}
import org.apache.avro.generic.GenericDatumWriter
import org.apache.avro.specific.SpecificDatumReader

// import hadoop stuff
import org.apache.hadoop.mapreduce.Job
import org.apache.hadoop.conf.Configuration

/**
 * Sequence Graph companion object for important IO methods. TODO: move this to
 * some other utility thing.
 */
object SequenceGraph {

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
    protected def writeRDDToParquet[RecordType <: IndexedRecord](
        things: RDD[RecordType], directory: String)
        (implicit m: Manifest[RecordType]) = {
        
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
        
        println("Writing out an RDD of %d %s".format(things.count, 
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
            
        println("RDD written")
                    
    }
    
    /**
     * Utility method to read a Parquet Avro format RDD. Handles the problem
     * that comes up when records are read in as the base specific record rather
     * than the actual code-generated type, which occurs when parquet-mr and
     * your code-generated types are in different jars.
     */
    def readRDDFromParquet[RecordType <: IndexedRecord](sc: SparkContext, 
        directory: String)
        (implicit m: Manifest[RecordType]): RDD[RecordType] = {
        
        // Every time we switch Avro schemas, we need a new Job. So we create
        // our own.
        val job = new Job(sc.hadoopConfiguration)
        
        // Configure Parquet to read the correct Avro thing. This may not
        // actually work to produce the right record type, but it's worth a
        // shot.
        ParquetInputFormat.setReadSupportClass(job, 
            classOf[AvroReadSupport[RecordType]])
            
        // Grab the configuration from the job, so things can load classes.
        val config = ContextUtil.getConfiguration(job)
        
        // Get the class of the record type. The type checker isn't too sure
        // this is going to work, but I'm pretty sure this is how manifests
        // work.
        val recordClass: java.lang.Class[RecordType] = m.erasure.asInstanceOf[
            java.lang.Class[RecordType]]
            
        // Read in the records as an RDD. We treat this as an RDD of the base
        // IndexedRecord type, which both the correct code-generated class and
        // the incorrect base specific record are instances of. We need to
        // specify the type for the method since Scala can't figure it out.
        val rdd: RDD[IndexedRecord] = sc.newAPIHadoopFile[java.lang.Void,
            RecordType,parquet.hadoop.ParquetInputFormat[RecordType]](directory, 
            classOf[ParquetInputFormat[RecordType]], classOf[Void], 
            recordClass, config)
            .map(p => p._2)
            
            
        // TODO: Check if we need this slow serialization loop to fix types.
        rdd.mapPartitions { (loadedIterator: Iterator[IndexedRecord]) =>
            // Do per-partition setup so we don't reflect in the inner loop.
            
            // Get the schema Field object for the requested record type   
            val recordSchemaField = recordClass.getDeclaredField("SCHEMA$")
            // Get the static value of this field. Argument is ignored for
            // static fields, so pass null. See <http://docs.oracle.com/javase/7
            // /docs/api/java/lang/reflect/Field.html#get%28java.lang.Object%29>
            val recordSchema = recordSchemaField.get(null).asInstanceOf[Schema]
            
            // Make a reader that produces things of the right type
            val datumReader = new SpecificDatumReader[RecordType](recordClass)
        
            loadedIterator.map { (loaded: IndexedRecord) =>
            
                // Serailize and de-serialize with Avro to get into the
                // appropriate type. Avro deserialization will load with the
                // thread context classloader (or the classloader of the class
                // we pass it), even though ParquetAvro deserialization won't.
                
                // This seems like a hack, but it's much easier to implement
                // than walking and fixing up the whole object graph from
                // ParquetAvro ourselves.
                
                // TODO: Fix the underlying bug and make ParquetAvro use the
                // thread context classloader.
               
                
                // Make an output stream that saves in a byte array
                val byteOutput = new ByteArrayOutputStream
                
                // Get an Avro Encoder to write to it, using the schema it was
                // read with.
                // TODO: Switch to binary encoding/decoding.
                val encoder = EncoderFactory.get.jsonEncoder(loaded.getSchema,
                    byteOutput)
                
                // Get a GenericDatumWriter to write to the encoder, using the
                // schema it was read with.
                val writer = new GenericDatumWriter[IndexedRecord](
                    loaded.getSchema)
                
                // Write
                writer.write(loaded, encoder)            
                encoder.flush()
                
                // Grab the bytes
                val bytes = byteOutput.toByteArray
                
                // Make an input stream
                val byteInput = new ByteArrayInputStream(bytes)
                
                // Make an Avro Decoder for the correct class, using our version of
                // the schema.
                val decoder = DecoderFactory.defaultFactory
                    .jsonDecoder(recordSchema, byteInput)
                
                try {
                    // Read
                    datumReader.read(null.asInstanceOf[RecordType], decoder)
                } catch {
                    case _: java.io.EOFException =>
                        // This will happen if you're trying to read with the
                        // wrong schema. TODO: fix schema migration.
                        throw new Exception("Check schema versions. EOF in: %s"
                            .format(new String(bytes)))
                }
            
            }
        
        }.cache
        
    }
}

/**
 * Represents a Sequence Graph (or component thereof) as a GraphX graph. The
 * nodes in the graph carry Sides, and the edges in the graph carry HasEdge
 * objects, which in turn carry AlleleGroups, Adjacencies, and Anchors.
 *    
 * Supports input/output to Parquet Avro, or dumping to any SequenceGraphWriter
 * in serial, and also some graph queries.
 */
class SequenceGraph(graph: Graph[Side, HasEdge]) {
    
    // Import fake static methods
    import SequenceGraph._
    
    /**
     * Extract all the Sides from this SequenceGraph.
     */
    def sides: RDD[Side] = {
        graph.vertices.map(_._2)
    }
    
    /**
     * Extract all the AlleleGroups from this SequenceGraph.
     */
    def alleleGroups: RDD[AlleleGroup] = {
        graph.edges.map(_.attr).flatMap {
            // Every AlleleGroup gets passed
            case AlleleGroupEdge(alleleGroup) => Some(alleleGroup)
            // Everything else gets thrown out
            case _ => None
        }
    }
    
    /**
     * Extract all the Adjacencies from this SequenceGraph.
     */
    def adjacencies: RDD[Adjacency] = {
        graph.edges.map(_.attr).flatMap {
            case AdjacencyEdge(adjacency) => Some(adjacency)
            case _ => None
        }
    }
    
    /**
     * Extract all the Anchors from this SequenceGraph.
     */
    def anchors: RDD[Anchor] = {
        graph.edges.map(_.attr).flatMap {
            case AnchorEdge(anchor) => Some(anchor)
            case _ => None
        }
    }
    
    /**
     * Make a SequenceGraph from separate RDDs of all the components.
     */
    def this(sidesRDD: RDD[Side], alleleGroupsRDD: RDD[AlleleGroup], 
        adjacenciesRDD: RDD[Adjacency], anchorsRDD: RDD[Anchor]) {
        
        this({
            // Make an RDD for Sides by ID
            val nodes: RDD[(Long, Side)] = sidesRDD.keyBy(_.id)
            
            // Make an RDD of all the different edge types wrapped up in HasEdge
            // objects. We have to map the implicit converters ourselves since
            // we don't have a magic RDD implicit converter mapper.
            val edges: RDD[HasEdge] = alleleGroupsRDD.map(AlleleGroup2HasEdge _) 
                .union(adjacenciesRDD.map(Adjacency2HasEdge _))
                .union(anchorsRDD.map(Anchor2HasEdge _))
        
            // Make an RDD of all the edges as GraphX Edge objects (not
            // SequenceGraph Edge objects) wrapping the HasEdge wrappers that
            // hold the actual Anchors, Adjacencies, and AlleleGroups.
            val graphEdges = edges.map { (edge) =>
                new org.apache.spark.graphx.Edge(edge.edge.left,
                    edge.edge.right, edge)
            }
            
            // Make and return the graph formed by these two RDDs.
            Graph(nodes, graphEdges)
        })
    }
    
    /**
     * Constructor to collect a bunch of SequenceGraphChunks together into a
     * SequenceGraph.
     */
    def this(parts: RDD[SequenceGraphChunk]) {
        // Pull out all the collections of parts, flatmap them into RDDs, and
        // use those RDDs.
        this(parts.flatMap(_.sides), parts.flatMap(_.alleleGroups), 
            parts.flatMap(_.adjacencies), parts.flatMap(_.anchors))
    }
    
    /**
     * Constructor to create a SequenceGraph in a given SparkContext from a
     * directory (or HDFS path) to which a SequenceGraph has previously been
     * saved.
     */
    def this(sc: SparkContext, directory: String) {
        // For some reason we need to SequenceGraph. this stuff even though we
        // imported it. And we can't even have an import before "this" here.
        this(SequenceGraph.readRDDFromParquet(sc, directory + "/Sides"), 
            SequenceGraph.readRDDFromParquet(sc, directory + "/AlleleGroups"), 
            SequenceGraph.readRDDFromParquet(sc, directory + "/Adjacencies"), 
            SequenceGraph.readRDDFromParquet(sc, directory + "/Anchors"))
    }
    
    /**
     * Combine this SequenceGraph with the other one, returning a SequenceGraph
     * with all the parts from both.
     */
    def union(other: SequenceGraph) = {
        // Make and return a new SequenceGraph holding all the unioned RDDs.
        new SequenceGraph(sides ++ other.sides, 
            alleleGroups ++ other.alleleGroups, 
            adjacencies ++ other.adjacencies, 
            anchors ++ other.anchors)
    }
    
    /**
     * Get the subgraph formed only from edges that overlap the given range. An
     * edge overlaps the range if it has at least one endpoint with a lower
     * bound position in the range, or if both its endpoints' lower bound
     * positions are on the same contig as the range and on opposite sides of
     * the range.
     *
     * Only vertices that have an edge overlapping the range will be preserved.
     *
     */
    def inRange(range: BaseRange): SequenceGraph = {
        new SequenceGraph(graph.subgraph(epred = {
            (triplet) =>
            // Pull out the endpoint position lower bounds. For things placed on
            // actual "reference" contigs, these will be real positions. For
            // things like novel inserts these will be on the contig the insert
            // went in to. TODO: Formalize this more.
            val leftPos = triplet.srcAttr.lowerBound
            val rightPos = triplet.dstAttr.lowerBound
            
            // Is the left Side in the range?
            val leftInRange = range contains leftPos 
            // Is the right Side in the range?
            val rightInRange = range contains rightPos
            // Are both Sides on the same contig as the range, but on opposite
            // sides of the range?
            val coversRange = range.between(leftPos, rightPos)
                
            // If any of those are true, we want this edge
            leftInRange || rightInRange || coversRange
        }).filter(preprocess = { (subgraph) =>
            // Label the subgraph vertices with their degrees. For some reason
            // we don't have a plain joinVertices that produces a graph with a
            // new type of vertex annotation, so we need to use this one and
            // pass a function.
            subgraph.outerJoinVertices(subgraph.degrees) {
                // We're replacing the vertex labels with the degree counts.
                // Things with degree 0 just don't appear, so we fill in for
                // them.
                (id, original, degree) => degree.getOrElse(0)
            }
        }, vpred = { (id: VertexId, degree: Int) =>
            // Only pass vertices with edges on them. We get the vertex
            // degrees by consulting the preprocessed graph, but actually
            // act on the original graph.
            degree > 0
        }))            
    }
    
    
    /**
     * An operator that unions this sequence graph with the other one and
     * returns the result.
     */
    def ++(other: SequenceGraph) = union(other)
    
    /**
     * Write this SequenceGraph to a series of Parquet directories.
     */
    def writeToParquet(directory: String) = {
    
        // Write each type of object in turn.
        writeRDDToParquet(sides, directory + "/Sides")
        writeRDDToParquet(adjacencies, directory + "/Adjacencies")
        writeRDDToParquet(alleleGroups, directory + "/AlleleGroups")
        writeRDDToParquet(anchors, directory + "/Anchors")
    
    }
    
    /**
     * Write this SequenceGraph to a (serial) SequenceGraphWriter. Collects
     * everything into driver node memory.
     *
     * TODO: Find a way to iterate over RDD items, fetching them lazily from
     * where they live.
     */
    def write(writer: SequenceGraphWriter) = {
    
        // Try doing an un-thread-safe thing over all the items.
        sides.collect.foreach(writer.writeSide _)
        adjacencies.collect.foreach(writer.writeAdjacency _)
        alleleGroups.collect.foreach(writer.writeAlleleGroup _)
        anchors.collect.foreach(writer.writeAnchor _)
        
    }
}

/**
 * Represents a bunch of Sequence Graph parts, such as might be returned from a
 * function that generates a small piece of a Sequence Graph. Can be stored in
 * an RDD.
 *
 * Immutable, so all modification methods return modified versions.
 */
class SequenceGraphChunk(sidesSeq: Seq[Side] = Nil,
    alleleGroupsSeq: Seq[AlleleGroup] = Nil,
    adjacenciesSeq: Seq[Adjacency] = Nil,
    anchorsSeq: Seq[Anchor] = Nil) {
    
    // We keep around all the sequence graph parts. Constructor arguments don't
    // magically become fields we can reference on other instances.
    val sides = sidesSeq
    val alleleGroups = alleleGroupsSeq
    val adjacencies = adjacenciesSeq
    val anchors = anchorsSeq
    
    
    
    /**
     * Combine this SequenceGraphChunk with the other one, returning a
     * SequenceGraphChunk with all the parts from both.
     */
    def union(other: SequenceGraphChunk) = {
        // Make and return a new SequenceGraphChunk holding all the unioned RDDs.
        new SequenceGraphChunk(sides ++ other.sides, 
            alleleGroups ++ other.alleleGroups, 
            adjacencies ++ other.adjacencies, 
            anchors ++ other.anchors)
    }
    
    /**
     * An operator that unions this sequence graph with the other one and
     * returns the result.
     */
    def ++(other: SequenceGraphChunk) = union(other)
    
    /**
     * Add some Sides to this SequenceGraphChunk.
     */
    def addSides(moreSides: Seq[Side]) = {
        new SequenceGraphChunk(sides ++ moreSides, alleleGroups, adjacencies, 
            anchors)
    }
        
    /**
     * Add some AlleleGroups to this SequenceGraphChunk.
     */
    def addAlleleGroups(moreAlleleGroups: Seq[AlleleGroup]) = {
        new SequenceGraphChunk(sides, alleleGroups ++ moreAlleleGroups, 
        adjacencies, anchors)
    }
    
    /**
     * Add some Adjacencies to this SequenceGraphChunk.
     */    
    def addAdjacencies(moreAdjacencies: Seq[Adjacency]) = {
        new SequenceGraphChunk(sides, alleleGroups, 
            adjacencies ++ moreAdjacencies, anchors)
    }
        
    /**
     * Add some Anchors to this SequenceGraphChunk.
     */
    def addAnchors(moreAnchors: Seq[Anchor]) = {
        new SequenceGraphChunk(sides, alleleGroups, adjacencies,
            anchors ++ moreAnchors)
    }
    
    /**
     * Add a single Side to a SequenceGraphChunk.
     */
    def +(side: Side) = addSides(List(side))
    
    /**
     * Add a single AlleleGroup to a SequenceGraphChunk.
     */
    def +(alleleGroup: AlleleGroup) = addAlleleGroups(List(alleleGroup))
    
    /**
     * Add a single Adjacency to a SequenceGraphChunk.
     */
    def +(adjacency: Adjacency) = addAdjacencies(List(adjacency))
    
    /**
     * Add a single Anchor to a SequenceGraphChunk.
     */
    def +(anchor: Anchor) = addAnchors(List(anchor))
    
    /**
     * Get the Side with the given ID, or None if no Side with that ID exists.
     * 
     * Just does a simple scan through all the Sides, so don't let your
     * SequenceGraphChunks become too big.
     */
    def getSide(id: Long) : Option[Side] = {
        sides.find(_.id == id)
    }
    
    /**
     * Get the `count` minimal-position-and-face Sides in this chunk for each
     * contig, in sorted order (smallest first).
     *
     * TODO: Implement proper k-select
     */
    def getMinimalSides(count: Int = 1) : Map[String, Seq[Side]] = {
        import scala.util.Sorting
    
        // Grab count lowest within each contig
        sides.groupBy(_.position.contig).mapValues { (value) => 
            Sorting.stableSort(value)
                .slice(0, count)
        }
        
    }
    
    /**
     * Get the `count` maximal-position-and-face Sides in this chunk for each
     * contig, in sorted order (largest first).
     *
     * TODO: Implement proper k-select
     */
    def getMaximalSides(count: Int = 1) : Map[String, Seq[Side]] = {
        import scala.util.Sorting
    
        // Grab count highest within each contig
        sides.groupBy(_.position.contig).mapValues { (value) => 
            Sorting.stableSort(value)
                .reverse
                .slice(0, count)
        }
        
    }
}

 
