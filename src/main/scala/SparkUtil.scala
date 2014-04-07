package edu.ucsc.genome

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import org.apache.spark.rdd.RDD

import org.apache.spark.graphx
import org.apache.spark.graphx._

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
    
        // How many partitions should we have? Graph misbehaves if the node and
        // edge RDDs have unequal numbers of partitions. To avoid a shuffle, we
        // coalesce down to the minimum number of partitions.
        val partitions = Math.min(nodes.partitions.size, edges.partitions.size)
        
        // Coalesce to an equal number of partitions, and make and return
        // the graph formed by these two RDDs.
        Graph(nodes.coalesce(partitions), edges.coalesce(partitions))
    }

}
