package edu.ucsc.genome

import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.avro.generic.IndexedRecord

// Import parquet
import parquet.hadoop.{ParquetOutputFormat, ParquetInputFormat}
import parquet.avro.{AvroParquetOutputFormat, AvroWriteSupport, 
                     AvroReadSupport, AvroParquetWriter}
import parquet.hadoop.util.{ContextUtil, ConfigurationUtil}

// We need to make Paths for Parquet output.
import org.apache.hadoop.fs.Path

// And we use a hack to get at the (static) schemas of generic things
import org.apache.avro.Schema

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
        val job = new Job()
        
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
        
        // Now correct the types
        rdd.map { 
            case record: RecordType =>
                // We got the correct type after all. Assume all its members are
                // correctly typed. No need to do anything.
                record
            case record: IndexedRecord =>
                // We got an incorrect record type, since Parquet refused to
                // load the right class. We need to correct the type not only of
                // this record but of all its IndexedRecord-extending members.
                fixRecordClass(config, record)
        }
            
    }
    
    /**
     * Given an object tree of IndexedRecords, assume that the root record is
     * actually the given record type, and correct it and all of its members.
     */
    def fixRecordClass[RecordType <: IndexedRecord](config: Configuration,
        record: IndexedRecord)(implicit m: Manifest[RecordType]): RecordType = {
        
        // Get this record's schema, which helpfully came with it
        val schema = record.getSchema
        
        // What class should it be? Guaranteed to be some class since it's an
        // IndexedRecord.
        val properClassName = schema.getFullName
        
        // Resolve that class in our whole classpath. Use the parquet-mr
        // ConfigurationUtil class that it should have used to look it up
        // originally. This will complain if you're demanding that this method
        // return one thing, but the thing you've deserialized is not one of
        // those.
        val properClass = ConfigurationUtil.getClassFromConfig(config, 
            properClassName, m.erasure)
            
        // Make a new instance of that, assuming that the class has a public no-
        // argument constructor, which the Avro code generator provides. Since
        // ConfigurationUtil already did runtime type checking, we can just
        // asInstanceOf here.
        val properInstance: RecordType = properClass.newInstance
            .asInstanceOf[RecordType]
        
        // How many fields are in this record's schema? Guaranteed to have a
        // fields slist since it's an IndexedRecord.
        val numFields = schema.getFields.size
        
        // Up-convert all the fields if necessary
        (0 until numFields).map { (index) =>
            (index, record.get(index) match {
                // Figure out what to do with each field's value
                case innerRecord: IndexedRecord =>
                    // If it's an IndexedRecord, make sure it has the correct
                    // actual type
                    fixRecordClass(config, innerRecord)
                case other: java.lang.Object =>
                    // Just pass through anything else
                    other
                // TODO: Fix up stuff inside of arrays, or enums.
            })
        }.foreach {
            case (index, value) =>
                // Set that field of the object to the fixed-up value.
                properInstance.put(index, value)
        }
        
        
        // Return the fixed-up IndexedRecord
        properInstance
        
        // For each field
            // Get the value
            // If the value is an IndexedRecord
                // Recurse on it, returning another IndexedRecord that's 
                // internally fixed.
            // Else just keep that value.
        // Look up and load the type we are supposed to be with getFullName
        // Make a new instance of it.
        // Populate each field with the value we found for that field.
        // Return the fixed up version.
        
    }

}

/**
 * Represents a Sequence Graph (or component thereof) as a series of Spark RDDs.
 * Supports input/output to Parquet Avro, or dumping to any SequenceGraphWriter
 * in serial, and also some graph queries.
 */
class SequenceGraph(sidesRDD: RDD[Side], alleleGroupsRDD: RDD[AlleleGroup], 
    adjacenciesRDD: RDD[Adjacency], anchorsRDD: RDD[Anchor]) {
    
    import SequenceGraph._
    
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

 
