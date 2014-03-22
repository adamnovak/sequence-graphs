package edu.ucsc.genome

/**
 * An implicitly-convertable-to class for Sequence Graph objects that have (or
 * are represented by) edges in the graph: Sites, Breakpoints, Abstractions,
 * AlleleGroups, Adjacencies, and Anchors.
 *
 * We do this instead of a type class since you can't have a list of type class
 * instances of different underlying types and polymorphically dispatch added
 * methods.
 *
 * Allows the edge member to be accessed in all of them without proper
 * (untagged) union types existing.
 *
 * Implicit conversions to case subclasses of this type exist in package.scala.
 */
trait HasEdge {
    def edge: Edge
}

/**
 * A class to say that AlleleGroups are edges.
 */
case class AlleleGroupEdge(owner: AlleleGroup) extends HasEdge {
    def edge = owner.edge
}

/**
 * A class to say that Adjacencies are edges.
 */
case class AdjacencyEdge(owner: Adjacency) extends HasEdge {
    def edge = owner.edge
}

/**
 * A class to say that Anchors are edges.
 */
case class AnchorEdge(owner: Anchor) extends HasEdge {
    def edge = owner.edge
}

/**
 * A class to say that Sites are edges.
 */
case class SiteEdge(owner: Site) extends HasEdge {
    def edge = owner.edge
}

/**
 * A class to say that Breakpoints are edges.
 */
case class BreakpointEdge(owner: Breakpoint) extends HasEdge {
    def edge = owner.edge
}

/**
 * A class to say that Generalizations are edges.
 */
case class GeneralizationEdge(owner: Generalization) extends HasEdge {
    def edge = owner.edge
}
