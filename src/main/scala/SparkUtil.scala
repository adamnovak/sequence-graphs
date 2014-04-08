package edu.ucsc.genome

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import org.apache.spark.rdd.RDD

import org.apache.spark.graphx
import org.apache.spark.graphx._

import java.io.{ByteArrayInputStream, ByteArrayOutputStream}

// Import parquet
import parquet.hadoop.{ParquetOutputFormat, ParquetInputFormat}
import parquet.avro.{AvroParquetOutputFormat, AvroWriteSupport, AvroReadSupport,
    AvroParquetWriter}
import parquet.hadoop.util.{ContextUtil, ConfigurationUtil}

// We need to make Paths for Parquet output.
import org.apache.hadoop.fs.Path

// And we use a hack to get at the (static) schemas of generic things
import org.apache.avro.Schema

import org.apache.avro.generic.IndexedRecord

// And more hacks to get things into the right type when Parquet Avro ignores
// our classpath.
import org.apache.avro.io.{EncoderFactory, DecoderFactory}
import org.apache.avro.generic.GenericDatumWriter
import org.apache.avro.specific.SpecificDatumReader

// import hadoop stuff
import org.apache.hadoop.mapreduce.Job
import org.apache.hadoop.conf.Configuration

import scala.reflect.ClassTag

/**
 * Object of utility methods for Spark/GraphX, extending it with operations it
 * does not natively support.
 *
 * Most operations require some amount of collecting things to the master.
 */
object SparkUtil {

    /**
     * Annotate every item in the given RDD with its index.
     */
    def zipWithIndex[T](rdd: RDD[T]): RDD[(T, Long)] = {
        
        // Count the number of items in each partition. Note that partition
        // iterators aren't infinite so we can use length on them.
        val partitionCounts = rdd.mapPartitions(i => Some(i.length).iterator)
        
        // Broadcast that
        val broadcastCounts = rdd.context.broadcast(partitionCounts.collect)
        
        rdd.mapPartitionsWithIndex({ (partition, iterator) =>
            // For each partition, count up the number of things before it to
            // get a base index.
            val baseIndex = broadcastCounts.value.take(partition).sum
        
            // Put base index + partition index with each item.
            iterator.zipWithIndex.map {
                case (value, index) => (value, index + baseIndex)
            }
        }, preservesPartitioning = true)
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
     * Make a graph from an RDD of vertices (keyed by ID) and an RDD of Edges.
     * You would think that you could just use GraphX's built-in Graph.apply to
     * do this (as it has the same signature), but that method has an
     * undocumented requirement that the vertex and edge RDDs have the same
     * number of partitions. This function allows vertex and edge RDDs with any
     * number of partitions, and coalesces to the minimum number.
     */
    def graph[VD: ClassTag, ED: ClassTag](nodes: RDD[(VertexId, VD)], 
        edges: RDD[org.apache.spark.graphx.Edge[ED]]): Graph[VD, ED] = {
    
        // How many partitions are in use in our input?
        val nodesBefore = nodes.partitions.size
        val edgesBefore = edges.partitions.size
    
        // How many partitions should we have? Graph misbehaves if the node and
        // edge RDDs have unequal numbers of partitions. To avoid a shuffle, we
        // coalesce down to the minimum number of partitions.
        val partitions = Math.min(nodes.partitions.size, edges.partitions.size)
        
        // Coalesce
        val nodesCoalesced = nodes.coalesce(partitions)
        val edgesCoalesced = edges.coalesce(partitions)
        
        // How many partitions are in use after?
        val nodesAfter = nodesCoalesced.partitions.size
        val edgesAfter = edgesCoalesced.partitions.size
        
        println("Coalesced to from %d/%d to %d/%d".format(nodesBefore, 
            edgesBefore, nodesAfter, edgesAfter))
        
        // Coalesce to an equal number of partitions, and make and return
        // the graph formed by these two RDDs.
        Graph(nodesCoalesced, edgesCoalesced)
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
    def writeRDDToParquet[RecordType <: IndexedRecord](
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
