package edu.ucsc.genome

import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.avro.generic.IndexedRecord

// Import parquet
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
 * Represents a Sequence Graph (or component thereof) as a series of Spark RDDs.
 */
class SequenceGraph(sidesRDD: RDD[Side], alleleGroupsRDD: RDD[AlleleGroup], 
    adjacenciesRDD: RDD[Adjacency], anchorsRDD: RDD[Anchor]) {
    
    // We keep around all the sequence graph parts. Constructor arguments don't
    // magically become fields we can reference on other instances.
    val sides = sidesRDD
    val alleleGroups = alleleGroupsRDD
    val adjacencies = adjacenciesRDD
    val anchors = anchorsRDD
    
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
     * An operator that unions this sequence graph with the other one and
     * returns the result.
     */
    def ++(other: SequenceGraph) = union(other)
    
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
        
        println("Writing out an RDD of %s".format(recordClass.getCanonicalName))
        
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

 
