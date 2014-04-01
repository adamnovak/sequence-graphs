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
     * Given a lower-level graph, create Sides and Edges for the next level up
     * according to this merging scheme.
     */
    def createNewLevel(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        (RDD[Side], RDD[HasEdge])
}

/**
 * Represents a not-actually-merged merging scheme. Just passes everything
 * through, with new IDs.
 */
case class Unmerged extends MergingScheme {
    def createNewLevel(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        (RDD[Side], RDD[HasEdge]) = {
        
        // IDs in the old graph must be distinct from IDs in this new graph.
        
        // Grab the SparkContext
        val sc = lowerLevel.edges.sparkContext
        
        // How many vertices are there?
        val vertexCount = lowerLevel.vertices.count
        
        // Get a new ID for each
        val sideIDStart: Long = ids.ids(vertexCount)
        
        // Get a new ID for the generalization from each.
        val generalizationIDStart = ids.ids(vertexCount)
        
        // Get a new ID for each existing edge.
        val edgeIDStart = ids.ids(lowerLevel.edges.count)
        
        // Get an ID to name the contig after. TODO: use level height.
        val contigID = ids.id
        
        // Assign each to a vertex.
        val annotations = SparkUtil.zipWithIndex(lowerLevel.vertices).map {
            case ((vertexID, vertex), index) => 
                (vertexID, index + sideIDStart)
        }
        
        // Join those annotations into the graph. TODO: find a way to do these
        // both as VertexRDDs.
        val annotatedGraph = lowerLevel.outerJoinVertices(annotations) {
            // Each vertex is now the vertex and the annotation, which we know
            // exists for every vertex.
            (vertexID, vertex, annotation) => (vertex, annotation.get)
        }
        
        val newSides = annotatedGraph.vertices.map { 
            // Fix up and return Side copies for each vertex
            case (vertexID, (side, newID)) =>
                // Make a copy of the Side and set its ID.
                val newSide = Side.newBuilder(side).setId(newID).build
                
                // Fix its Position to be on a contig for this new level.
                // TODO: do this in order somehow, for compression of runs.
                newSide.position.contig = "merged%d".format(contigID)
                newSide.position.base = newID
                
                // Return it
                newSide
        }
        
        
        val newEdges = SparkUtil.zipWithIndex(annotatedGraph.triplets).map { 
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

        // Add Generalizations from old Sides to new Sides.
        val generalizations: RDD[HasEdge] = SparkUtil.zipWithIndex(
            annotatedGraph.vertices).map { 
            
            // Make a new Generalization from the old vertex to the new one,
            // using the appropriate edge ID based on the index of the vertex.
            case ((vertexID, (side, newID)), index) =>
                new Generalization(new Edge(vertexID, newID, 
                    index + generalizationIDStart))
        }
        
        // Send back the Sides and the edges (both ones in the upper level and
        // ones leading into it).
        (newSides, newEdges.union(generalizations))
    }
}

/**
 * Represents the context-driven merging scheme. TODO: Not going to implement
 * now.
 */
case class ContextDriven extends MergingScheme {
    def createNewLevel(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        (RDD[Side], RDD[HasEdge]) = {
        // TODO: Implement this
        val sc = lowerLevel.edges.sparkContext
        (sc.parallelize(Nil), sc.parallelize(Nil))
    }
}

/**
 * Represents the symmetric merging scheme requiring the given amount of context
 * on both sides.
 *
 * You have to iterate to a fixed point.
 */
case class Symmetric(context: Int) extends MergingScheme {
    def createNewLevel(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        (RDD[Side], RDD[HasEdge]) = {
        // TODO: Implement this
        val sc = lowerLevel.edges.sparkContext
        (sc.parallelize(Nil), sc.parallelize(Nil))
    }
}

/**
 * Represents the symmetric merging scheme requiring the given amount of context
 * on one side.
 */
case class NonSymmetric(context: Int) extends MergingScheme {
    def createNewLevel(lowerLevel: Graph[Side, HasEdge], ids: IDSource): 
        (RDD[Side], RDD[HasEdge]) = {
        // TODO: Implement this
        val sc = lowerLevel.edges.sparkContext
        (sc.parallelize(Nil), sc.parallelize(Nil))
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
            // next iteration.
            levelGraph = makeGraph(newSides, newEdges)
            
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
                        levelNumber) => true
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
