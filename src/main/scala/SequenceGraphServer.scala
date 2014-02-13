package edu.ucsc.genome

import scala.collection.JavaConversions._

/**
 * A class that implements the Avro-defined SequenceGraphAPI, allowing a
 * SequenceGraph to be queried via RPC calls.
 */
class SequenceGraphServer(graph: SequenceGraph) extends SequenceGraphAPI {

    /**
     * Get a subgraph of our SequenceGraph as a paginated response. Callable via
     * RPC.
     */
    def search(contig: String, start: Long, end: Long): GraphElementResponse = {
    
        // Subset the graph
        val subgraph = graph.inRange(new BaseRange(contig, start, end))
        
        // Make an RDD of the parts
        val elements = subgraph.sides.map(new GraphElement(_))
            .union(subgraph.alleleGroups.map(new GraphElement(_)))
            .union(subgraph.adjacencies.map(new GraphElement(_)))
            .union(subgraph.anchors.map(new GraphElement(_)))
            
        // Send the entire thing at once.
        new GraphElementResponse(elements.collect.toList, elements.count, "")
    }
    
    /**
     * Get the next page of a paginated response. Callable via RPC.
     */
    def getNextPage(token: String): GraphElementResponse = {
        
        // We always put everything on the first page, so there will never be
        // any more pages.
        new GraphElementResponse(List[GraphElement](), 0, "")
    }
    
}
