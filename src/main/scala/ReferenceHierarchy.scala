package edu.ucsc.genome

import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.graphx._

import org.apache.commons.io.FileUtils
import java.io._

import com.esotericsoftware.kryo._
import com.esotericsoftware.kryo.io.{Input, Output}
import com.twitter.chill.ScalaKryoInstantiator

import scala.reflect._

/**
 * Represents the kind of reference structure that a level in a reference
 * hierarchy is (i.e. what merging scheme is used to produce it).
 */
trait MergingScheme {
    /**
     * Given a lower-level graph and a place to get new unique IDs, create Sides
     * and Edges for the next level up according to this merging scheme.
     */
    def createNewLevel(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        (RDD[Side], RDD[HasEdge]) = {
        
        // Call scheme-specific logic to decide what nodes to merge.
        val annotatedGraph = annotate(lowerLevel, ids)

        // Make Sides for each pair of sets of merged Sides.
        // Collect Sides by new ID, and sort by (contig, base)
        val setsToMerge = annotatedGraph.vertices.groupBy {
            // Group by newID
            case (vertexId, (side, newID)) => newID
        }.map {
            case (newID, annotatedSides) => 
                (newID, annotatedSides.map {
                    // Turn into (newID, list of sides)
                    case (vertexId, (side, _)) => side
                }
                // And make sure that list of sides is sorted by (contig, base)
                .sortBy(_.position.contig)
                .sortBy(_.position.base))
        }
        
        // Get an ID to name the contig after on the new layer. TODO: use level
        // height.
        val contigID = ids.id
        
        // Create all the merged Sides and with their IDs
        val newSides = setsToMerge.map { 
            case (newID, sides) =>
                
                // Make the new Side, adopting the Face of the first thing in
                // the list. There must be one Side at least in each list, or
                // the list wouldn't exist.
                val firstSide = sides(0)
                
                // Make a copy of the Side and set its ID.
                val newSide = Side.newBuilder(firstSide).setId(newID).build
                
                // Fix its Position to be on a contig for this new level. Just
                // rename its contig and use its original base number, since
                // we're guaranteed the partner can do the same.
                newSide.position.contig = "%s-merged%d"
                    .format(newSide.position.contig, contigID)
                // Keep its Face as whatever it had originally. The other side
                // of that base will also be the first Side of its group sorted
                // by (contig, base) if we are always merging both sides of
                // bases, so we will have a partner that is our opposite face.
                
                // TODO: Maybe re-number bases here to have something to do with
                // Site edge IDs?
                
                // TODO: fix up lowerBounds and other Side fields.
                
                // Return it
                newSide
        }
        
        // Find one Site edge to represent every new merged pair of Sides, and
        // fix it up. We assume that both Sides of each Site edge are going to
        // be merged into corrsponding Sides of a new Site edge.
        val siteTriplets = annotatedGraph.triplets.filter { triplet => 
            // Select only the Sites
            triplet.attr.isInstanceOf[SiteEdge]
        }
        // Group with the other Sites that are getting merged into this new Site.
        .groupBy { triplet =>
            Math.min(triplet.srcAttr._2, triplet.dstAttr._2)
        }
        // Take an arbitrary triplet to represent each Site. It will have its ID
        // and endpoints fixed up.
        .map { case (_, triplets) => triplets(0) }
        
        // Find one Breakpoint edge to represent all the breakpoints between a
        // pair of new Sides, in any order.
        val breakpointTriplets = annotatedGraph.triplets.filter { triplet => 
            // Select only the Breakpoints
            triplet.attr.isInstanceOf[BreakpointEdge]
        }
        // Group with the other Breakpoints between these new Sites.
        .groupBy { triplet =>
            // Pull out the two new IDs we are connecting
            val srcID = triplet.srcAttr._2
            val dstID = triplet.dstAttr._2
            
            // Put into a tuple in order.
            (Math.min(srcID, dstID), Math.max(srcID, dstID))
        }
        // Take an arbitrary triplet to represent each connection. It will have
        // its ID and endpoints fixed up.
        .map { case (_, triplets) => triplets(0) }
        
        // Put the filtered Site and Breakpoint triplets together.
        val keptEdges = siteTriplets.union(breakpointTriplets)
            
        // Get a new ID for each kept edge.
        val edgeIDStart = ids.ids(keptEdges.count)
            
        // Copy all the Site and Breakpoint triplets and change their IDs and
        // endpoints.
        val newEdges = SparkUtil.zipWithIndex(keptEdges).map {
            case (triplet, index) =>
                // Fix up and return edge copies for each triplet
                val clone = triplet.attr.clone
                
                // Set the left node ID to the new ID for the left node
                clone.edge.left = triplet.srcAttr._2
                // And similarly for the right node
                clone.edge.right = triplet.dstAttr._2
                
                // And give it a new ID
                clone.edge.id = index + edgeIDStart
                
                // Return it
                clone
        }
        
        // Get a new ID for the generalization from each vertex on the lower
        // level.
        val generalizationIDStart = ids.ids(lowerLevel.vertices.count)

        // Add Generalizations from old Sides to new Sides.
        val generalizations: RDD[HasEdge] = SparkUtil.zipWithIndex(
            annotatedGraph.vertices).map { 
            
            // Make a new Generalization from the old vertex to the new one,
            // using the appropriate edge ID based on the index of the vertex.
            case ((vertexId, (side, newID)), index) =>
                new Generalization(new Edge(index + generalizationIDStart, 
                    vertexId, newID))
        }
        
        // Send back the Sides and the edges (both ones in the upper level and
        // ones leading into it).
        (newSides, newEdges.union(generalizations))
    
    }
        
    /**
     * Given a lower-level graph and a place to get new unique IDs, create an
     * annotated graph where each Side is paired with a number identifying the
     * new Side it ought to be merged into. These numbers are not unique; sides
     * sharing a number wil be merged together.
     */
    def annotate(lowerLevel: Graph[Side, HasEdge], ids: IDSource):
        Graph[(Side, Long), HasEdge]
}

/**
 * Represents a not-actually-merged merging scheme. Just passes everything
 * through, with new IDs.
 */
case class Unmerged extends MergingScheme {
    def annotate(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        Graph[(Side, Long), HasEdge] = {
        
        // IDs in the old graph must be distinct from IDs in this new graph.

        // How many new IDs do we need?
        val vertexCount = lowerLevel.vertices.count 

        // Get a new ID for each vertex
        val sideIDStart: Long = ids.ids(vertexCount)
        
        // Attach indices to vertices
        val indexed = SparkUtil.zipWithIndex(lowerLevel.vertices)
        
        // Assign each to a vertex.
        val annotations = indexed.map {
            case ((vertexId, vertex), index) => 
                (vertexId, index + sideIDStart)
        }
        
        // Join those annotations into the graph. TODO: find a way to do these
        // both as VertexRDDs.
        lowerLevel.outerJoinVertices(annotations) {
            // Each vertex is now the vertex and the annotation, which we know
            // exists for every vertex.
            (vertexId, vertex, annotation) => (vertex, annotation.get)
        }
        
        
    }
}

/**
 * Represents the context-driven merging scheme. TODO: Not going to implement
 * now.
 */
case class ContextDriven extends MergingScheme {
    def annotate(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        Graph[(Side, Long), HasEdge] = {
        
        // TODO: Implement this
        null
    }
}

/**
 * Represents the symmetric merging scheme requiring the given amount of context
 * on both sides.
 *
 * You have to iterate to a fixed point.
 */
case class Symmetric(context: Int) extends MergingScheme {
    def annotate(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        Graph[(Side, Long), HasEdge] = {
        
        // TODO: Implement this
        null
    }
}

/**
 * Represents the nonsymmetric merging scheme requiring the given amount of
 * context on one side.
 */
case class NonSymmetric(context: Int) extends MergingScheme {
    def annotate(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        Graph[(Side, Long), HasEdge] = {
        
        println(("Annotating lower level with %d vertex partitions and %d " + 
            "edge partitions").format(lowerLevel.vertices.partitions.size, 
            lowerLevel.edges.partitions.size))
        
        println("Annotating")
        
        // Go tag each Side with all the contexts of exactly the right length
        // that it appears in.
        val contextGraph = runSearch(lowerLevel, context)
        
        println("Grouping by context")
        
        // Group Sides by context string. We assume that no more than a
        // reasonable number of Sides share a context string.
        val idsByContext = contextGraph.vertices.flatMap { 
            case (id, (side, contexts)) => contexts.map((_, id))
        }.groupByKey
        
        // Make a new graph where each Side is in a connected component with all
        // other Sides sharing a context with it. The sides are the sides from
        // the original graph, but the edges are made from the grouped ID sets.
        
        println("Making context edges")
        
        // Put together the edges
        val contextEdges = idsByContext.flatMap { case (context, sideIds) =>
            // We can't get here with an empty ID list. We only want to add
            // edges between the first side that has a context and any
            // subsequent side sharing that context..
            val first = sideIds(0)
            val rest = sideIds.drop(1)
            
            rest.map { other =>
                // Make these edges as true, since they connect things in the
                // same component.
                new org.apache.spark.graphx.Edge(first, other, true)
            }
        }
        
        println("Making complementary edges")
        
        // Add all the Sites as edges, marking connected components that must be
        // complementary.
        val complementEdges = lowerLevel.edges.flatMap { edge =>
            // Look at the edge attribute type.
            edge.attr match {
                case siteEdge: SiteEdge => 
                    // We need an edge with attribute false for this Site, so
                    // its two Sides end up in complementary components.
                    Some(new org.apache.spark.graphx.Edge(edge.srcId,
                        edge.dstId, false))
                case _ =>
                    // We need to drop this edge
                    None
            }
        }
        
        println("Making problem graph")
        
        // Put together the complementary connected components input graph.
        // Union the two edge RDDs and coalesce to the number of partitions in
        // the vertexRDD we are re-using. Should preserve the vertex index
        // structure.
        val problemGraph =  Graph(lowerLevel.vertices, 
            contextEdges.union(complementEdges)
            .coalesce(lowerLevel.vertices.partitions.size))
        
        
        println("Solving Complementary Connected Components")
        
        // Solve, so that each component has at most one complementary
        // component, and Sides connect complementary components. Now we have a
        // graph where the vertex IDs are side IDs, and the vertex attributes
        // are the components that those side IDs ought to be merged into. If
        // two nodes share the same attribute, they need to merge.
        val componentGraph = ComplementaryConnectedComponents.run(problemGraph)
        
        println("Finding free IDs")
        
        // What's the maximum vertex value in the component graph? We can't need
        // more than that many + 1 IDs. TODO: Assumes no negative IDs.
        val maxComponent = componentGraph.vertices.map(_._2).fold(0)(Math.max _)
        
        // Re-number the components so they don't conflict with the original
        // vertex id values. What base number should we use for them? We ened to
        // reserve enough so that when we offset the max component by this
        // amount, we'll be safely under the number we reserved.
        val componentIdStart = ids.ids(maxComponent + 1)
        
        println("Zipping")
        
        // This holds the annotated vertices for our final graph. Since we have
        // the same vertices we can use the cheap join.
        val annotatedVertices = lowerLevel.vertices
            .innerZipJoin(componentGraph.vertices) { (id, side, component) =>
                // Annotate each side with an unused ID corresponding to its
                // component
                (side, componentIdStart + component)            
            }
        
        // Stick the annotated vertices back together with their original edges.
        // We can use Graph here since we know we didn't change the partition
        // counn for annotatedVertices since we split it off of lowerLevel.
        val toReturn = Graph(annotatedVertices, lowerLevel.edges)
        
        println("Annotated")
        toReturn
    }
    
    /**
     * Run a Pregel search to find all the upstream contexts of length exactly
     * length (not counting the base that the node belongs to, which is also
     * included) for each node.
     */
    protected def runSearch(graph: Graph[Side, HasEdge], length: Int): 
        Graph[(Side, List[String]), HasEdge] = {
        
        // First, run the search to completion, so each node has a list of
        // SearchStates with no breadcrumbs left.
        
        println("Creating search graph")
        
        // Set up the initial graph.
        val searchGraph: Graph[(Side, List[SearchState]), HasEdge] = graph
            .mapVertices((id, side) => (side, Nil) )
        
        // What should we start out with for each node? A search state
        // programmed to cross to the other side of its Site, then go down depth
        // and back up depth again.
        val initialMessage = List(new SearchState(length * 2 + 1))
        
        println("Running Pregel")
        
        // Run Pregel for enough iterations for every search state to get to the
        // other side of its Site, then down to depth and back up again.
        val pregelGraph = Pregel(searchGraph, initialMessage, 
            length * 2 + 1)(vertexProgram _, sendMessage _, messageCombiner _)
        
        println("Preparing answers")
        
        // Map so each node has a list of Strings, on the strand corresponding
        // to upstream.
        pregelGraph.mapVertices { (id, attr) => 
            // Join all the characters each state found into a string. Keep the
            // Side as is.
            (attr._1, attr._2.map(_.characters.mkString))
        }
    }
    
    /**
     * Figure out what to do at each vertex in the context search.
     */
    protected def vertexProgram(id: VertexId, attr: (Side, List[SearchState]),
        msgSum: List[SearchState]): (Side, List[SearchState]) = {
        
        // The search states on this node are replaced by the ones we just got.
        (attr._1, msgSum)
    }
    
    
    /**
     * Figure out what messages should be sent along each edge in the context
     * search. Runs once per edge per iteration.
     */
    protected def sendMessage(
        edge: EdgeTriplet[(Side, List[SearchState]), HasEdge]): 
        Iterator[(VertexId, List[SearchState])] = {
     
        // We should send one searchstate in each direction for each search
        // state at a node that hasn't reached its end.
        
        // Get the messages from the source vertex, and tag with the
        // destination ID
        val fromSrc: List[(VertexId, List[SearchState])] = edge.srcAttr._2
            .map(_.getMessagesOver(edge.srcId, edge))
            .map((edge.dstId, _))
            
        // Get the messages from the destination vertex, and tag with the
        // source ID
        val fromDst: List[(VertexId, List[SearchState])] = edge.dstAttr._2
            .map(_.getMessagesOver(edge.dstId, edge))
            .map((edge.srcId, _))
        
        // Make an iterator for both of these lists.
        fromSrc.iterator ++ fromDst.iterator
    }
    
    /**
     * Combine two messages coming into a node in the context search.
     */
    def messageCombiner(a: List[SearchState], b: List[SearchState]): 
        List[SearchState] = {
        
        // Just concatenate the lists
        a ++ b
    }

    
}


/**
 * A modified version of the connected components algorithm where you have two
 * types of edges. One type of edge links nodes into the same component. The
 * other type of edge links "complementary" components. Components are merged
 * such that each component has at most one complementary component.
 *
 * In this particular problem, the first type of edges are the shared-context
 * edges, and the second type are Site edges. If two Sides of one Site share a
 * component, the other two Sides need to also share a component.
 *
 * Based on the GraphX library code. See <https://github.com/amplab/graphx/
 * blob/master/graphx/src/main/scala/org/apache/spark/graphx/lib/ConnectedCompon
 * ents.scala>.
 */
object ComplementaryConnectedComponents {

    // What type should messages be? Also the type of vertex attributes. Holds
    // component id and complementary component ID.
    type Message = (Long, Long)

    /**
     * Compute the "Complementary Connected Components" component membership of
     * each vertex.  Return a graph where each vertex has the lowest vertex ID
     * of any vertex in its component, and each component has at most one
     * complementary component.
     *
     * Edge attributes should be true for edges that link vertices in the same
     * component, and false for edges that link vertices that must be in
     * complementary components.
     *
     */
    def run[VD: ClassTag](graph: Graph[VD, Boolean]):
        Graph[VertexId, Boolean] = {
        
        // We implement this by having each node track the minimum vertex ID in
        // its component, and the minimum vertex ID in its complementary
        // component.
        
        // Label the graph vertices with their IDs as their min component value,
        // and the max long as their complementary component value.
        val startingGraph = graph.mapVertices { 
            case (id: VertexId, _) => (id, Long.MaxValue)
        }
        
        // Start with an initial message of the max possible values (which will
        // do nothing).
        val initialMessage = (Long.MaxValue, Long.MaxValue)
        
        // Run Pregel to tag each vertex with its component and complementary
        // component. Run until we stop sending any messages.
        val solvedGraph = Pregel(startingGraph, initialMessage, 
            activeDirection=EdgeDirection.Either)(vertexProgram, sendMessage, 
            messageCombiner)
          
        // Return the graph with just the actual component IDs, dropping the
        // complementary IDs.
        solvedGraph.mapVertices {
            case (id, attr) => attr._1
        }
    } 
    
    /**
     * Given a vertex ID, its vertex attribute, and the message it got, return
     * its new attribute.
     */
    protected def vertexProgram(id: VertexId, attr: Message,
        message: Message): Message = {
        
        // Just use the message combiner to get the min of everything.
        messageCombiner(attr, message)
    }
    
    /**
     * Given an edge triplet, return an iterator of the destination IDs and
     * message values to send over the edge.
     */
    protected def sendMessage(edge: EdgeTriplet[Message, Boolean]): 
        Iterator[(VertexId, Message)] = {
        
        // Grab the source and destination tuples.
        val src = edge.srcAttr
        val dst = edge.dstAttr
        
        // Make a List of (id, message) pairs to send.
        var toSend: List[(VertexId, Message)] = Nil
        
        if(edge.attr) {
            // This edge links nodes in the same component.
            
            if(src._1 > dst._1 || src._2 > dst._2) {
                // We need to tell the source about the destination.
                toSend ::= ((edge.srcId, dst))
            }
            
            if(dst._1 > src._1 || dst._2 > src._2) {
                // We need to tell the destination about the source.
                toSend ::= ((edge.dstId, src))
            }
        } else {
            // This edge links nodes in complementary components.
            
            if(src._1 > dst._2 || src._2 > dst._1) {
                // We need to tell the source about the destination, backwards.
                toSend ::= ((edge.srcId, (dst._2, dst._1)))
            }
            
            if(dst._1 > src._2 || dst._2 > src._1) {
                // We need to tell the destination about the source, backwards.
                toSend ::= ((edge.dstId, (src._2, src._1)))
            }
        }
        
        // Send all the messages
        toSend.iterator
    }
    
    /**
     * Given two messages (each of which is a component ID and a complementary
     * component ID), return a message with the minimum component ID and the
     * minimum complementary component ID.
     */
    def messageCombiner(a: Message, b: Message): Message = {
        // Just (manually) map min over the zipped tuples.
        (Math.min(a._1, b._1), Math.min(a._2, b._2))
    }
}

/**
 * Represents the state of a search in a reference structure. Wants to traverse
 * edges alternating between Breakpoints and Sites, splitting up at every branch
 * point, going out to a specified depth of Sites, and collecting the base
 * sequences from the paths traversed.
 *
 * Each search traverses the Site it starts at, and ends at the *opposite* side,
 * with the contexts of the *opposite* side. Since each Site has two Sides, each
 * Side still ends up with a list of search states that traversed its upstream
 * context.
 *
 * breadcrumbs is a stack of the path we took, in terms of node IDs we need to
 * go back to.
 *
 * characters is a stack of string characters we encountered.
 *
 * depthRemaining is the total number of edges left to traverse. After
 * traversing the initial Site edge, this is 2 * the number of characters we
 * still have to get, and is even if we have to take a Breakpoint and odd if we
 * have to take a Site.
 *
 */
class SearchState(val depthRemaining: Int, val breadcrumbs: List[Long] = Nil, 
    val characters: List[Char] = Nil) extends Serializable  {
    
    
    
    /**
     * Get all the messages (SearchStates, should only ever be 1) that should be
     * sent form the given vertex over the given edge.
     */
    def getMessagesOver(sender: VertexId, 
        edge: EdgeTriplet[(Side, List[SearchState]), HasEdge]): 
        List[SearchState] = {
        
        if(depthRemaining > 0) {
            // Traverse edges downwards.
            edge.attr match {
                case siteEdge: SiteEdge =>
                    if(depthRemaining % 2 == 1) {
                        // Use Site edges when odd
                        
                        // Work out what character to add: edge's or reverse
                        // complement? We take the forward character when we
                        // traverse the edge *backwards*, since we're builing
                        // our string from end to start.
                        val charToAdd = sender match {
                            case id if id == edge.srcId =>
                                // Pull it out as a char first. We're going
                                // forwards, so take the reverse complement.
                                siteEdge.owner.base(0).reverseComplement
                            case _ =>
                                // We're going backwards, so leave the character
                                // alone after we pull it out.
                                siteEdge.owner.base(0)
                        }
                        
                        // Leave a breadcrumb so we come back here, if
                        // applicable.
                        val newBreadcrumbs = characters match {
                            // Don't leave a breadcrumb on the very first node.
                            // So we'll return to the opposite Side of our Site.
                            case Nil => breadcrumbs
                            // Leave breadcrumbs when crossing all subsequent
                            // Sites.
                            case _ => sender :: breadcrumbs
                        }
                        
                        // Add in this character and continue down.
                        List(new SearchState(depthRemaining - 1, newBreadcrumbs,
                            charToAdd :: characters))
                    } else {
                        // Don't use Site edges when even
                        Nil
                    }
                case breakpointEdge : BreakpointEdge =>
                    if(depthRemaining % 2 == 0) {
                        // Use breakpoint edges when even.
                    
                        // Send a copy along this edge, reducing
                        // depthRemaining and adding a breadcrumb.
                        List(new SearchState(depthRemaining - 1,
                            sender :: breadcrumbs, characters))
                    } else {
                        Nil
                    }
            }
        } else {
            // We're going back up.
            breadcrumbs match {
                case head :: rest if head == edge.otherVertexId(sender) =>
                    // This is the edge we came down. Go back up it and pop it
                    // off the stack.
                    List(new SearchState(depthRemaining, rest, characters))
                case _ =>
                    // Don't send any other messages.
                    Nil
            }
            
        }
        
    }
    
}

/**
 * Represents a reference hierarchy composed of reference structures at
 * different levels. Defined by a bottom-level index, a merging scheme for
 * building each level, and a possibly finished, or possibly null, graph of
 * Sides, Sites, Adjacencies, and Generalizations (with no negative vertex IDs).
 */
class ReferenceHierarchy(sc: SparkContext, var index: FMDIndex) {
    
    // Make an IDSource for generating sequential IDs. TODO: Fix up for
    // load/save.
    val source = new IDSource(0)
    
    // Keep a Seq of ReferenceStructure levels, from bottom to top
    var levels: List[ReferenceStructure] = Nil
    
    // Keep around a level-labeled version of the graph.
    var labeledGraph: Graph[(Side, Int), HasEdge] = null
    
    // Keep around the graph without level labels
    var graph: Graph[Side, HasEdge] = null
    
    /**
     * Load a ReferenceHierarchy from the given path, using the given
     * SparkContext. The saved hierarchy must contain at least one level.
     */
    def this(sc: SparkContext, path: String) = {
        // Set up with a null index.
        this(sc, null.asInstanceOf[FMDIndex])
        
        // Load our graph.
        
        // Load the Sides
        val sides: RDD[Side] = SparkUtil
            .readRDDFromParquet(sc, path + "/Sides")
        // And all the edges
        val siteEdges: RDD[Site] = SparkUtil
            .readRDDFromParquet(sc, path + "/Sites")
        val breakpointEdges: RDD[Breakpoint] = SparkUtil
            .readRDDFromParquet(sc, path + "/Breakpoints")
        val generalizationEdges: RDD[Generalization] = SparkUtil
            .readRDDFromParquet(sc, path + "/Generalizations")
            
        // Key the Sides by ID
        val keyedSides = sides.keyBy(_.id)
            
        // Make the edges into HasEdges
        val haveEdges: RDD[HasEdge] = siteEdges.map(Site2HasEdge)
            .union(breakpointEdges.map(Breakpoint2HasEdge))
            .union(generalizationEdges.map(Generalization2HasEdge))
            
        // And then into Graphx Edges
        val edges = haveEdges.map { hasEdge =>
            new org.apache.spark.graphx.Edge(hasEdge.edge.left, 
                hasEdge.edge.right, hasEdge)
        }
        
        // And make the graph
        graph = SparkUtil.graph(keyedSides, edges)
        
        // Deserialize levels with Kryo. TODO: Ensure registrator stuff has
        // happened.
        val kryo = (new ScalaKryoInstantiator).newKryo
        val input = new Input(new FileInputStream(path + "/levels.dat"))
        
        // Load our labeled levels.
        levels = kryo.readObject(input, classOf[List[ReferenceStructure]])
        
        input.close
        
        // Grab the index from the first one, which we assume exists.
        index = levels.head.getIndex
        
        // TODO: Try to do this without mutable state, and without needing to do
        // loads of stuff in an argument to this().
    }
    
    
    /**
     * Save this ReferenceHierarchy to the specified directory.
     */
    def save(path: String) = {
        // Make a File 
        val directory = new File(path)
        
        if(directory.exists) {
            if(directory.isDirectory) {
                // Delete the directory if it exists.
                FileUtils.deleteDirectory(directory)
            } else {
                // Delete it if it's a file, too.
                directory.delete
            }
        }
        
        // Make it again
        directory.mkdir
        
        // Save with Kryo
        val kryo = (new ScalaKryoInstantiator).newKryo
        val output = new Output(new FileOutputStream(path + "/levels.dat"))
        kryo.writeObject(output, levels);
        output.close
        
        // Save the unlabeled graph RDDs. We will unfortunately have to do
        // without the labels when we reload, since we can only save RDDs of
        // Avro stuff.
        
        // Grab the vertices as Sides
        val sides: RDD[Side] = labeledGraph.vertices.map(_._2._1)
        // And the edges
        val siteEdges: RDD[Site] = labeledGraph.edges
            .map(_.attr)
            .flatMap {
                case edge: SiteEdge => Some(edge.owner)
                case _ => None
            }
        val breakpointEdges: RDD[Breakpoint] = labeledGraph.edges
            .map(_.attr)
            .flatMap {
                case edge: BreakpointEdge => Some(edge.owner)
                case _ => None
            }
        val generalizationEdges: RDD[Generalization] = labeledGraph.edges
            .map(_.attr)
            .flatMap {
                case edge: GeneralizationEdge => Some(edge.owner)
                case _ => None
            }
        
        // Write everything    
        SparkUtil.writeRDDToParquet(sides, path + "/Sides")
        SparkUtil.writeRDDToParquet(siteEdges, path + "/Sites")
        SparkUtil.writeRDDToParquet(breakpointEdges, path + "/Breakpoints")
        SparkUtil.writeRDDToParquet(generalizationEdges, path +
            "/Generalizations")       
        
    }
    
    /**
     * Build a completely new graph by merging, starting from our bottom-level
     * contigs. Also build the levels list.
     */
    def initialize(schemes: Seq[MergingScheme]) = {
        
        // Make an empty RDD of sides for our graph, paired with the hierarchy
        // levels they live at, starting at 0 on the bottom.
        var sides: RDD[(Side, Int)] = sc.parallelize(Nil)
        // And another empty RDD of various types of edges
        var edges: RDD[HasEdge] = sc.parallelize(Nil)
        
        for(contig <- index.contigs) {
            // Make sides and edges for each contig
            val (newSides, newEdges) = makeContigGraph(contig)
            
            // Add them to our growing collections of graph parts.
            sides = sides.union(newSides.map((_, 0)))
            edges = edges.union(newEdges)
        }
        
        // Keep a graph of just each level by itself, without level labels,
        // starting from this bottom level.
        var levelGraph = makeGraph(sides.map(_._1), edges)
        
        for((scheme, schemeIndex) <- schemes.zipWithIndex) {
            println("Merging level %d...".format(schemeIndex + 1))
            
            // Make new sides and edges from merging on the last graph. TODO:
            // this will eventually need to be able to map things, for some
            // schemes.
            val (newSides, newEdges) = scheme.createNewLevel(levelGraph, source)
            
            // Make a new graph with just the new stuff, and keep it for the
            // next iteration. But make sure we don't include ghost nodes at the
            // lower ends of the generalization edges. TODO: Don't mix and
            // immediately unmix these things.
            levelGraph = makeGraph(newSides, newEdges).subgraph(vpred = {
                (vertexId, vertex) => vertex != null
            })
            
            println("Made level with %d sides and %d edges".format(
                newSides.count, newEdges.count))
            
            // Roll these sides and edges into the collection we are going to
            // use to create the big unioned graph. We can't just union our
            // small graphs. Make sure to tag the sides with their level,
            // starting at 1 for the first merged level.
            sides = sides.union(newSides.map((_, schemeIndex + 1)))
            edges = edges.union(newEdges)
        }
        
        // Now we've created levels for all of our schemes. Make our overall
        // graph, with level labels. We need to explain that the Sides to get
        // the IDs from are the first things in the tuples. We also need to
        // cache the graph, so we don't re-do all that work.
        labeledGraph = makeGraph(sides, edges,
            {(tuple: (Side, Int)) => tuple._1}).cache()
        
        // Make the final graph with the labels stripped
        graph = labeledGraph.mapVertices {
            (id, data) => data._1
        }
        
        
        
        // Now we made the graph; we need to look at it and make the actual
        // levels we use for mapping.
        
        // Initialize the levels list with just the bottom-level structure
        levels = List(new StringReferenceStructure(index))
        
        for(levelNumber <- 1 until schemes.size + 1) {
            // For each reference structure we need to build on top, from the
            // bottom up...
            
            println("Making level %d...".format(levelNumber))
            
            // Make a new ReferenceStructure on top of the previous one. The
            // contig passed here should never get used since we don't let the
            // ReferenceStructure auto-generate any Positions; they're all
            // specified by the merged graphs.
            val newStructure = new CollapsedReferenceStructure(levels.head, 
                "level%d".format(levelNumber))
            
            // Find all the tripples for Generalization edges feeding into that
            // level.
            val generalizations = labeledGraph.subgraph { triplet =>
                (triplet.attr, triplet.dstAttr._1.position.face, 
                    triplet.dstAttr._2) match {
                    // Break out based on edge type, destination face, and
                    // destination layer.
                    
                    // We found a generalization going into the appropriate
                    // level, on the left side of something.
                    case (generalization: GeneralizationEdge, Face.LEFT, 
                        nodeLevel) => nodeLevel == levelNumber
                    // We found something else
                    case _ => false
                }
            }
            
            // Group Positions to merge by destination Position. TODO: do we
            // need to use longs here instead of Positions as keys?
            val positionsToMerge = generalizations.triplets.map { triplet =>
                (triplet.dstAttr._1.position, triplet.srcAttr._1.position)
            }.groupByKey
        
            // Put the new structure on top of the stack.
            levels = newStructure :: levels
        }
        
        // Now we've built all the reference structure levels. Flip them over,
        // so they go from bottom to top.
        levels = levels.reverse
        
    }
    
    /**
     * Given the names of new contigs added to the bottom level, re-build the
     * higher levels, and replace graph with an update graph.
     * 
     * Graph may actually be empty, and we are building it initially.
     */
    def update(newContigs: Seq[String]) = {
    
        // TODO: Implement after the working prototype.
    
        // Build a graph for the new contigs at the bottom level.
        
        // Union it in with our existing graph.
        
        // For each level:
        
            // Fix up BWT intervals so we can map to/search the lower level.
            
                // We need this for identifying nodes to merge with.
                
            // Add to-merge edges connecting corresponding Sides of lower-level
            // nodes to merge.
            
            // Find the connected components of the lower-level nodes along
            // these edges.
            
            // Find all the upper-level nodes connected to each component by
            // Generalization edges.
            
            // Take the first pair of Positions from the upper-level nodes as
            // the Positions to use for IDs for each connected component. Use
            // the remainder as aliases. If no Positions were found, create new
            // ones.
            
            // Make a new set of upper-level nodes, one per merged connected
            // component.
            
            // Add in the Adjacency edges between the merged nodes.
            
            // Add in the Generalization edges up from what was merged.
            
            // Add in Generalization edges up to the next higher level, for
            // finding IDs when we merge things down here.
            
            // Throw out all the to-merge edges, and the old upper-level nodes.
            // We don't need to add any to resolve the case of multiple
            // generalization destinations because that will be resolved
            // automatically when we replace the next- higher level.
            
        // When we've gotten to the top, we are done.
    }
    
    /**
     * Make the vertices and edges for the bottom-level graph for a given
     * contig.
     */
    def makeContigGraph(contig: String): (RDD[Side], RDD[HasEdge]) = {
        // Keep a growing list of Sides
        var sides: List[Side] = Nil
        // And a growing list of Sites/Breakpoints
        var edges: List[HasEdge] = Nil
        
        for(base <- 1L until index.contigLength(contig) + 1) {
            // Make positions for the base, accounting for 1-based indexing
            val leftPosition = new Position(contig, base, Face.LEFT)
            val rightPosition = new Position(contig, base, Face.RIGHT)
            
            // Make Sides for the base
            val leftSide = new Side(source.id, leftPosition, false,
                leftPosition)
            val rightSide = new Side(source.id, rightPosition, false,
                rightPosition)
            
            // Make the Site edge with the base for this position that we pull
            // from the FMDIndex. TODO: it might be much easier to just display
            // the whole contig and then iterate through that.
            edges = new Site(new Edge(source.id, leftSide.id, rightSide.id), 
                index.display(leftPosition).toString) :: edges
                
            if(sides != Nil) {
                // And the Breakpoint edge to the previous Site, if any
                edges = new Breakpoint(new Edge(source.id, sides.head.id, 
                    leftSide.id), false) :: edges
            }
            
            // Put our Sides in the list (right side last so it can be grabbed
            // on the next iteration).
            sides = rightSide :: leftSide :: sides
            
        }
        
        // Turn the Lists we built into RDDs and return them.
        (sc.parallelize(sides), sc.parallelize(edges))
    }
    
    /**
     * Perform the necessary gymnastics to make a GraphX graph out of some sides
     * and some edges. Sides may be annotated by being placed in tuples or
     * something, so you need to provide a function to unpack them.
     */
    def makeGraph[T](sides: RDD[T], edges: RDD[HasEdge],
        sideFunction: T => Side)(implicit t: ClassTag[T]): Graph[T, HasEdge] = {
        
        // Make an RDD for Sides by ID
        val nodes: RDD[(Long, T)] = sides.keyBy(side => sideFunction(side).id)
        
        // Make an RDD of all the edges as GraphX Edge objects (not
        // SequenceGraph Edge objects) wrapping the HasEdge wrappers.
        val graphEdges = edges.map { (edge) =>
            new org.apache.spark.graphx.Edge(edge.edge.left,
                edge.edge.right, edge)
        }
        
        // Coalesce to an equal number of partitions, and make and return
        // the graph formed by these two RDDs.
        val toReturn = SparkUtil.graph(nodes, graphEdges)
        
        println("Made %d/%d graph".format(toReturn.vertices.partitions.size, 
            toReturn.edges.partitions.size))
        
        toReturn
        
    }
    
    /**
     * Make a graph from unannotated Sides and their edges.
     */
    def makeGraph(sides: RDD[Side], edges: RDD[HasEdge]):
        Graph[Side, HasEdge] = {
        
        makeGraph(sides, edges, identity)
    }
    
}

/**
 * Class to dump a hierarchy graph to a file in GraphViz format.
 */
class GraphvizWriter(file: String) {
    // We want to write files.
    import java.io._
    
    // We need to write escaped strings. See
    // <http://stackoverflow.com/a/9914380/402891> and
    // <http://commons.apache.org/proper/commons-
    // lang/apidocs/org/apache/commons/lang3/StringEscapeUtils.html>
    import org.apache.commons.lang.StringEscapeUtils.escapeJava
    
    // Open the file for writing. See
    // <http://www.tutorialspoint.com/scala/scala_file_io.htm>
    val graphWriter = new java.io.PrintWriter(new File(file))
    
    // Write a header
    graphWriter.write("digraph hierarchy {\n")
    
    // Set up edge styles
    graphWriter.write("edge [arrowsize=\"0\"];\n")
    
    // And node styles
    graphWriter.write("node [shape=\"point\"];\n")
    
    /**
     * Write out a graph.
     */
    def writeGraph(graph: Graph[Side, HasEdge]) = {
        // Write all the Sides
        graph.vertices.map(_._2).collect.foreach(writeSide _)
        // Write all the edges by type
        graph.edges.map(_.attr).collect.foreach {
            case edge: SiteEdge => writeSite(edge)
            case edge: BreakpointEdge => writeBreakpoint(edge)
            case edge: GeneralizationEdge => writeGeneralization(edge)
            case _ => throw new Exception("Unhandled edge type")
        }
        
        close()
    }
    
    /**
     * Write a graph where each vertex is labeled with a subgraph.
     */
    def writeSubgraphs(graph: Graph[(Side, Int), HasEdge]) = {
        
        // Put all the Sides into seqs, keyed by subgraph, and collect.
        val verticesBySubgraph = graph.vertices.map { 
            case (id, (side, subgraph)) =>
                (subgraph, side)
        }.groupByKey.collect
        
        verticesBySubgraph.foreach { case (subgraph, vertices) =>
            // Start a subgraph. Make sure to start it's name with "cluster".
            // See <http://www.graphviz.org/Gallery/directed/cluster.html>
            graphWriter.write("subgraph cluster%d {\n".format(subgraph))
            // Style it so we can see it.
            graphWriter.write("style=filled;\ncolor=lightgrey;\n")
            // Label it with the level
            graphWriter.write("label = \"Level %d\"\n".format(subgraph + 1))
            
            // Do each vertex
            vertices.foreach(writeSide _)
            
            // End the subgraph
            graphWriter.write("}\n")
        }
        
        // Write all the edges by type
        graph.edges.map(_.attr).collect.foreach {
            case edge: SiteEdge => writeSite(edge)
            case edge: BreakpointEdge => writeBreakpoint(edge)
            case edge: GeneralizationEdge => writeGeneralization(edge)
            case _ => throw new Exception("Unhandled edge type")
        }
        
        close()
    }
    
    // Implementations for writing edges and nodes.
    
    def writeSite(edge: SiteEdge) : Unit = {
        // What should we label it?
        val label = edge.owner.base
        
        // Make an edge for every Site.
        // TODO: escape labels
        graphWriter.write(
            "%d -> %d [label=\"%s\",arrowsize=\"1\",color=\"#ff0000\"];\n"
            .format(edge.edge.left, edge.edge.right,
            escapeJava(label)))
    }
    
    def writeBreakpoint(edge: BreakpointEdge) : Unit = {
        // Make an edge for every Breakpoint
        graphWriter.write("%d -> %d;\n".format(
            edge.edge.left, edge.edge.right))
    }
    
    def writeGeneralization(edge: GeneralizationEdge) : Unit ={
        // Make an edge for every Generalization
        graphWriter.write("%d -> %d [arrowsize=\"1\",color=\"#00ff00\"];\n"
            .format(edge.edge.left, edge.edge.right))
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
            
        // Connect them with an edge
        graphWriter.write(
            "%d -> L%d [dir=\"none\",style=\"dotted\"];\n"
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
