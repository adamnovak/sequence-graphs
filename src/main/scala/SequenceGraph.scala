package edu.ucsc.genome

import java.io.{ByteArrayOutputStream, ByteArrayInputStream, DataOutputStream,
    DataInputStream}

import scala.collection.JavaConversions._

import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.graphx._

/**
 * Represents a Sequence Graph (or component thereof) as a GraphX graph. The
 * nodes in the graph carry Sides, and the edges in the graph carry HasEdge
 * objects, which in turn carry AlleleGroups, Adjacencies, and Anchors.
 *    
 * Supports input/output to Parquet Avro, or dumping to any SequenceGraphWriter
 * in serial, and also some graph queries.
 */
class SequenceGraph(graph: Graph[Side, HasEdge]) {
    
    /**
     * Extract all the Sides from this SequenceGraph.
     */
    def sides: RDD[Side] = {
        graph.vertices.map(_._2)
    }
    
    /**
     * Extract all the AlleleGroups from this SequenceGraph.
     */
    def alleleGroups: RDD[AlleleGroup] = {
        graph.edges.map(_.attr).flatMap {
            // Every AlleleGroup gets passed
            case AlleleGroupEdge(alleleGroup) => Some(alleleGroup)
            // Everything else gets thrown out
            case _ => None
        }
    }
    
    /**
     * Extract all the Adjacencies from this SequenceGraph.
     */
    def adjacencies: RDD[Adjacency] = {
        graph.edges.map(_.attr).flatMap {
            case AdjacencyEdge(adjacency) => Some(adjacency)
            case _ => None
        }
    }
    
    /**
     * Extract all the Anchors from this SequenceGraph.
     */
    def anchors: RDD[Anchor] = {
        graph.edges.map(_.attr).flatMap {
            case AnchorEdge(anchor) => Some(anchor)
            case _ => None
        }
    }
    
    /**
     * Make a SequenceGraph from separate RDDs of all the components.
     */
    def this(sidesRDD: RDD[Side], alleleGroupsRDD: RDD[AlleleGroup], 
        adjacenciesRDD: RDD[Adjacency], anchorsRDD: RDD[Anchor]) {
        
        this({
            // Make an RDD for Sides by ID
            val nodes: RDD[(Long, Side)] = sidesRDD.keyBy(_.id)
            
            // Make an RDD of all the different edge types wrapped up in HasEdge
            // objects. We have to map the implicit converters ourselves since
            // we don't have a magic RDD implicit converter mapper.
            val edges: RDD[HasEdge] = alleleGroupsRDD.map(AlleleGroup2HasEdge _) 
                .union(adjacenciesRDD.map(Adjacency2HasEdge _))
                .union(anchorsRDD.map(Anchor2HasEdge _))
                
            // Make an RDD of all the edges as GraphX Edge objects (not
            // SequenceGraph Edge objects) wrapping the HasEdge wrappers that
            // hold the actual Anchors, Adjacencies, and AlleleGroups.
            val graphEdges = edges.map { (edge) =>
                new org.apache.spark.graphx.Edge(edge.edge.left,
                    edge.edge.right, edge)
            }
            
            // How many partitions should we have? Graph misbehaves (when we do
            // inRange) if the node and edge RDDs have unequal numbers of
            // partitions. To avoid a shuffle, we coalesce down to the minimum
            // number of partitions.
            val partitions = Math.min(nodes.partitions.size,
                graphEdges.partitions.size)
            
            // Coalesce to an equal number of partitions, and make and return
            // the graph formed by these two RDDs.
            Graph(nodes.coalesce(partitions), graphEdges.coalesce(partitions))
        })
    }
    
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
        this(SparkUtil.readRDDFromParquet(sc, directory + "/Sides"), 
            SparkUtil.readRDDFromParquet(sc, directory + "/AlleleGroups"), 
            SparkUtil.readRDDFromParquet(sc, directory + "/Adjacencies"), 
            SparkUtil.readRDDFromParquet(sc, directory + "/Anchors"))
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
     * Get the subgraph formed only from edges that overlap the given range. An
     * edge overlaps the range if it has at least one endpoint with a lower
     * bound position in the range, or if both its endpoints' lower bound
     * positions are on the same contig as the range and on opposite sides of
     * the range.
     *
     * Only vertices that have an edge overlapping the range will be preserved.
     *
     */
    def inRange(range: BaseRange): SequenceGraph = {
        new SequenceGraph(graph.subgraph(epred = {
            (triplet) =>
            
            if(triplet == null || triplet.srcAttr == null || 
                triplet.dstAttr == null) {
                // Pre-emptively fail if we managed to get nulls somehow. TODO:
                // we shouldn't ever get nulls in here. ScalaTest is getting
                // them when running tests somehow.
                false
            } else {
            
                // Pull out the endpoint position lower bounds. For things
                // placed on actual "reference" contigs, these will be real
                // positions. For things like novel inserts these will be on the
                // contig the insert went in to. TODO: Formalize this more.
                val leftPos = triplet.srcAttr.lowerBound
                val rightPos = triplet.dstAttr.lowerBound
                
                // Is the left Side in the range?
                val leftInRange = range contains leftPos 
                // Is the right Side in the range?
                val rightInRange = range contains rightPos
                // Are both Sides on the same contig as the range, but on
                // opposite sides of the range?
                val coversRange = range.between(leftPos, rightPos)
                    
                // If any of those are true, we want this edge
                leftInRange || rightInRange || coversRange
            }
        }).filter(preprocess = { (subgraph) =>
            // Label the subgraph vertices with their degrees. For some reason
            // we don't have a plain joinVertices that produces a graph with a
            // new type of vertex annotation, so we need to use this one and
            // pass a function.
            subgraph.outerJoinVertices(subgraph.degrees) {
                // We're replacing the vertex labels with the degree counts.
                // Things with degree 0 just don't appear, so we fill in for
                // them.
                (id, original, degree) => degree.getOrElse(0)
            }
        }, vpred = { (id: VertexId, degree: Int) =>
            // Only pass vertices with edges on them. We get the vertex
            // degrees by consulting the preprocessed graph, but actually
            // act on the original graph.
            degree > 0
        }))            
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
        SparkUtil.writeRDDToParquet(sides, directory + "/Sides")
        SparkUtil.writeRDDToParquet(adjacencies, directory + "/Adjacencies")
        SparkUtil.writeRDDToParquet(alleleGroups, directory + "/AlleleGroups")
        SparkUtil.writeRDDToParquet(anchors, directory + "/Anchors")
    
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

 
