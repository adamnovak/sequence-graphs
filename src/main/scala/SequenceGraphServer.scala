package edu.ucsc.genome


import scala.reflect.ClassTag
import scala.collection.JavaConversions._
import org.apache.spark.rdd.RDD

/**
 * Represents an ongoing upload of a graph to a client. Keeps track of the RDD
 * being sent, the next page to send, and the cached count of the RDD.
 */
class OngoingTransfer(toSend: RDD[GraphElement], page: Int, total: Long) {
    
    /**
     * Start a new transfer of the given RDD.
     */
    def this(rdd: RDD[GraphElement]) = {
        this(rdd, 0, rdd.count)
    }

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
     * Count the number of pages remaining to send, after this one. When this
     * reaches 0, the transfer is finished.
     */
    def pagesRemaining: Int = toSend.partitions.size - (page + 1)
    
    /**
     * Check to see if this is the last page.
     */
    def isDone: Boolean = pagesRemaining == 0
    
    /**
     * Get the graph elements for this page.
     */
    def getPage: List[GraphElement] = {
        getPartition(toSend, page).collect.toList
    }
    
    /**
     * Make a GraphElementResponse to send this page to the client.
     * Automatically sets the nextPageToken to a unique value if this was not
     * the last page, or "" if this was the last page.
     */
    def response: GraphElementResponse = {
        new GraphElementResponse(getPage, total, 
            if(isDone) "" else java.util.UUID.randomUUID.toString)
    }
    
    /**
     * Make an OngoingTransfer to represent the need to send the next page. Make
     * sure to check isDone on it before using it.
     */
    def advance: OngoingTransfer = {
        new OngoingTransfer(toSend, page + 1, total)
    }
}

/**
 * A class that implements the Avro-defined SequenceGraphAPI, allowing a
 * SequenceGraph to be queried via RPC calls. Thread safe by virtue of only
 * letting one operation be happening at a time.
 */
class SequenceGraphServer(graph: SequenceGraph) extends SequenceGraphAPI {

    // Keep track of the state of every in-progress transfer of pages to
    // clients.
    var nextPages: Map[String, OngoingTransfer] = Map.empty

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
                
            // Make a transfer to send pages to the client.
            val transfer = new OngoingTransfer(elements)
            
            // Get the response for the first page
            val response = transfer.response
            
            if(!transfer.isDone) {
                // Remember that we have to continue the transfer of the next
                // page in response to the token.
                nextPages += ((response.nextPageToken, transfer.advance))
            }
                
            // Send the response
            response
        }
    }
    
    /**
     * Get the next page of a paginated response. Callable via RPC.
     */
    def getNextPage(token: String): GraphElementResponse = {
        this.synchronized {
            if(nextPages contains token) {
            
                // Grab the ongoing transfer state.
                val transfer = nextPages(token)    
                
                // Revoke the token
                nextPages -= token  
                
                // Work out what to send
                val response = transfer.response
                
                if(!transfer.isDone) {
                    // Remember that we have to continue the transfer of the
                    // next page in response to the token.
                    nextPages += ((response.nextPageToken, transfer.advance))
                }
                
                // Send the response
                response
            } else {
                // They're asking for something we ran out of.
                throw new BadRequestError
            }
        }
    }
    
}
