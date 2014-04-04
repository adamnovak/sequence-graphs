package edu.ucsc.genome

import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.graphx._


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
            case (vertexID, (side, newID)) => newID
        }.map {
            case (newID, annotatedSides) => 
                (newID, annotatedSides.map {
                    // Turn into (newID, list of sides)
                    case (vertexID, (side, _)) => side
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
                
                // Fix its Position to be on a contig for this new level.
                // TODO: do this in order somehow, for compression of runs.
                newSide.position.contig = "merged%d".format(contigID)
                newSide.position.base = newID
                // Keep its Face as whatever it had originally. The other side
                // of that base will also be the first Side of its group sorted
                // by (contig, base) if we are always merging both sides of
                // bases, so we will have a partner that is our opposite face.
                
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
            case ((vertexID, (side, newID)), index) =>
                new Generalization(new Edge(index + generalizationIDStart, 
                    vertexID, newID))
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

        // Get a new ID for each vertex
        val sideIDStart: Long = ids.ids(lowerLevel.vertices.count)
        
        // Assign each to a vertex.
        val annotations = SparkUtil.zipWithIndex(lowerLevel.vertices).map {
            case ((vertexID, vertex), index) => 
                (vertexID, index + sideIDStart)
        }
        
        // Join those annotations into the graph. TODO: find a way to do these
        // both as VertexRDDs.
        lowerLevel.outerJoinVertices(annotations) {
            // Each vertex is now the vertex and the annotation, which we know
            // exists for every vertex.
            (vertexID, vertex, annotation) => (vertex, annotation.get)
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
        
        // TODO: Implement this
        Unmerged().annotate(lowerLevel, ids)
    }
    
    /**
     * Run a Pregel search to find all the upstream contexts of length exactly
     * length (not counting the base that the node belongs to, which is also
     * included) for each node.
     */
    def runSearch(graph: Graph[Side, HasEdge], length: Int): 
        Graph[List[String], HasEdge] = {
        
        // First, run the search to completion, so each node has a list of
        // SearchStates with no breadcrumbs left.
        
        // Set up the initial graph.
        val searchGraph: Graph[List[SearchState], HasEdge] = graph
            .mapVertices((_, _) => Nil)
        
        // What should we start out with for each node? A search state
        // programmed to cross to the other side of its Site, then go down depth
        // and back up depth again.
        val initialMessage = List(new SearchState(length * 2 + 1))
        
        // Run Pregel for enough iterations for every search state to get to the
        // other side of its Site, then down to depth and back up again.
        val pregelGraph = Pregel(searchGraph, initialMessage, 
            length * 2 + 1)(vertexProgram _, sendMessage _, messageCombiner _)
        
        // Map so each node has a list of Strings, on the strand corresponding
        // to upstream.
        pregelGraph.mapVertices { (id, states) => 
            // Join all the characters each state found into a string.
            states.map(_.characters.mkString)
        }
    }
    
    /**
     * Figure out what to do at each vertex.
     */
    def vertexProgram(id: VertexId, attr: List[SearchState],
        msgSum: List[SearchState]): List[SearchState] = {
        
        // The search states on this node are the ones we just got.
        msgSum
    }
    
    
    /**
     * Figure out what messages should be sent along each edge. Runs once per
     * edge.
     */
    def sendMessage(edge: EdgeTriplet[List[SearchState], HasEdge]): 
        Iterator[(VertexId, List[SearchState])] = {
     
        // We should send one searchstate in each direction for each search
        // state at a node that hasn't reached its end.
        
        // Get the messages from the source vertex, and tag with the
        // destination ID
        val fromSrc: List[(VertexId, List[SearchState])] = edge.srcAttr
            .map(_.getMessagesOver(edge.srcId, edge))
            .map((edge.dstId, _))
            
        // Get the messages from the destination vertex, and tag with the
        // source ID
        val fromDst: List[(VertexId, List[SearchState])] = edge.dstAttr
            .map(_.getMessagesOver(edge.dstId, edge))
            .map((edge.srcId, _))
        
        // Make an iterator for both of these lists.
        fromSrc.iterator ++ fromDst.iterator
    }
    
    /**
     * Combine two messages coming into a node.
     */
    def messageCombiner(a: List[SearchState], b: List[SearchState]): 
        List[SearchState] = {
        
        // Just concatenate the lists
        a ++ b
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
    val characters: List[Char] = Nil)  {
    
    
    
    /**
     * Get all the messages (SearchStates, should only ever be 1) that should be
     * sent form the given vertex over the given edge.
     */
    def getMessagesOver(sender: VertexId, 
        edge: EdgeTriplet[List[SearchState], HasEdge]) : List[SearchState] = {
        
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
class ReferenceHierarchy(sc: SparkContext, index: FMDIndex,
    schemes: Seq[MergingScheme], var graph: Graph[Side, HasEdge] = null) {
    
    // Find the next available ID: the max in the graph we got, or 0.
    // Assumes the graph we got doesn't use any negative IDs.
    var nextAvailableID: Long = graph match {
        case null => 0
        case _ => graph.vertices.map(_._1)
            // Add in all the edge IDs too
            .union(graph.edges.map(_.attr.edge.id))
            // Get 1 after the largest, or 0.
            .fold(-1)((a, b) => if(a > b) a else b) + 1
    }
     
    // Make an IDSource for generating sequentila IDs. 
    val source = new IDSource(nextAvailableID)
    
    // Keep a Seq of ReferenceStructure levels, from bottom to top
    var levels: List[ReferenceStructure] = Nil
    
    /**
     * Build a completely new graph by merging, starting from our bottom-level
     * contigs. Also build the levels list.
     */
    def initialize = {
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
            // Make new sides and edges from merging on the last graph. TODO:
            // this will eventually need to be able to map things, for some
            // schemes.
            val (newSides, newEdges) = scheme.createNewLevel(levelGraph, source)
            
            // Make a new graph with just the new stuff, and keep it for the
            // next iteration. But make sure we don't include ghost nodes at the
            // lower ends of the generalization edges. TODO: Don't mix and
            // immediately unmix these things.
            levelGraph = makeGraph(newSides, newEdges).subgraph(vpred = {
                (vertexID, vertex) => vertex != null
            })
            
            // Roll these sides and edges into the collection we are going to
            // use to create the big unioned graph. We can't just union our
            // small graphs. Make sure to tag the sides with their level,
            // starting at 1 for the first merged level.
            sides = sides.union(newSides.map((_, schemeIndex + 1)))
            edges = edges.union(newEdges)
        }
        
        // Now we've created levels for all of our schemes. Make our overall
        // graph, with level labels. We need to explain that the Sides to get
        // the IDs from are the first things in the tuples.
        val labeledGraph = makeGraph(sides, edges,
            {(tuple: (Side, Int)) => tuple._1})
        
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
            
            println("=========================================================")
            println("MAKING LEVEL %d".format(levelNumber))
            
            // Make a new ReferenceStructure on top of the previous one. The
            // contig passed here should never get used since we don't let the
            // ReferenceStructure auto-generate any Positions; they're all
            // specified by the merged graphs.
            val newStructure = new CollapsedReferenceStructure(levels.head, 
                "level%d".format(levelNumber))
            
            // TODO: Finish this
            
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
        
            // Collect to the master
            val positionsCollected = positionsToMerge.collect
            
            
            for((destination, sources) <- positionsCollected) {
                println("Creating new position: %s".format(destination))
                sources.foreach(source => println("From %s".format(source)))
                
                // Add a new base with the correct destination Position
                // subsuming the source Positions of each group.
                newStructure.addBase(sources, destination)
            }
            
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
                index.display(leftPosition)) :: edges
                
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
        
        // How many partitions should we have? Graph misbehaves if the node and
        // edge RDDs have unequal numbers of partitions. To avoid a shuffle, we
        // coalesce down to the minimum number of partitions.
        val partitions = Math.min(nodes.partitions.size,
            graphEdges.partitions.size)
        
        // Coalesce to an equal number of partitions, and make and return
        // the graph formed by these two RDDs.
        Graph(nodes.coalesce(partitions), graphEdges.coalesce(partitions))
    }
    
    /**
     * Make a graph from unannotated Sides and their edges.
     */
    def makeGraph(sides: RDD[Side], edges: RDD[HasEdge]):
        Graph[Side, HasEdge] = {
        
        makeGraph(sides, edges, identity)
    }
    
}
