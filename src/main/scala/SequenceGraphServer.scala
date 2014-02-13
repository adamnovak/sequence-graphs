package edu.ucsc.genome


import scala.reflect.ClassTag
import scala.collection.JavaConversions._
import org.apache.spark.rdd.RDD

/**
 * A class that implements the Avro-defined SequenceGraphAPI, allowing a
 * SequenceGraph to be queried via RPC calls. Thread safe by virtue of only
 * letting one operation be happening at a time.
 */
class SequenceGraphServer(graph: SequenceGraph) extends SequenceGraphAPI {

    // Keep track of the next partition to send for every next page request.
    var nextPages: Map[String, (RDD[GraphElement], Int)] = Map.empty

    /**
     * Grab the specified partition of the given RDD as its own RDD.
     */
    def getPartition[T](rdd: RDD[T], partition: Int)
        (implicit tag: ClassTag[T]): RDD[T] = {
        
        rdd.mapPartitionsWithIndex { case (index: Int, data: Iterator[T]) =>
            if(index == partition) {
                // Keep everything from the correct partition.
                data
            } else {
                // Throw out everything in the other partitions.
                None.iterator
            }
        }
    }

    /**
     * Get a subgraph of our SequenceGraph as a paginated response. Callable via
     * RPC.
     */
    def search(contig: String, start: Long, end: Long): GraphElementResponse = {
        this.synchronized {
    
            // Subset the graph
            val subgraph = graph.inRange(new BaseRange(contig, start, end))
            
            // Make an RDD of the parts
            val elements = subgraph.sides.map(new GraphElement(_))
                .union(subgraph.alleleGroups.map(new GraphElement(_)))
                .union(subgraph.adjacencies.map(new GraphElement(_)))
                .union(subgraph.anchors.map(new GraphElement(_)))
                
            // Make a new string token
            val token: String = java.util.UUID.randomUUID.toString
            
            // Say the next page is partition 1 (if it exists)
            val entry = (token, (elements, 1))
            nextPages += entry
                
            // Send partition 0
            new GraphElementResponse(getPartition(elements, 0).collect.toList, 
                elements.count, token)
        }
    }
    
    /**
     * Get the next page of a paginated response. Callable via RPC.
     */
    def getNextPage(token: String): GraphElementResponse = {
        this.synchronized {
            if(nextPages contains token) {
            
                // See what we're supposed to be sending
                val (elements, partition) = nextPages(token)    
                
                // Revoke the token
                nextPages -= token  
                
                // Get the list to send back
                val list = getPartition(elements, partition).collect.toList
                
                // Make a new string token
                val nextToken: String = list match {
                    // If we ran out justs end them to the forever empty token.
                    case Nil => ""
                    case _ => {
                        // Generate a new token
                        val generated = java.util.UUID.randomUUID.toString
                        
                        // Say the next page is the next partition (if it
                        // exists)
                        val entry = (generated, (elements, partition + 1))
                        nextPages += entry
                        
                        generated
                    }
                }
                
                // Send requested partition
                new GraphElementResponse(list, elements.count, nextToken)
            } else {
                // They're asking for something we ran out of.
                // TODO: Don't suddenly say there were 0 results...
                new GraphElementResponse(List(), 0, "") 
            }
        }
    }
    
}
