package edu.ucsc.genome

import org.apache.spark.rdd.RDD

/**
 * Represents a Sequence Graph (or component thereof) as a series of Spark RDDs.
 */
class SequenceGraph(sidesRDD: RDD[Side], alleleGroupsRDD: RDD[AlleleGroup], 
    adjacenciesRDD: RDD[Adjacency], anchorsRDD: RDD[Anchor]) {
    
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
    
    
}

/**
 * Represents a bunch of Sequence Graph parts, such as might be returned from a
 * function that generates a small piece of a Sequence Graph. Can be stored in
 * an RDD.
 */
class SequenceGraphChunk(sidesSeq: Seq[Side], alleleGroupsSeq: Seq[AlleleGroup],
    adjacenciesSeq: Seq[Adjacency], anchorsSeq: Seq[Anchor]) {
    
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
    
}

/**
 * Provides an implicit conversion to extend RDDs of SequenceGraphChunks to
 * easily be turned into SequenceGraphs.
 */

 
