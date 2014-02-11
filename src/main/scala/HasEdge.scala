package edu.ucsc.genome

/**
 * An implicitly-cohnvertable-to class for Sequence Graph objects that have (or
 * are represented by) edges in the graph: AlleleGroups, Adjacencies, and
 * Anchors.
 *
 * Allows the edge member to be accessed in all of them without proper
 * (untagged) union types existing.
 */
trait HasEdge {
    def edge: Edge
}
