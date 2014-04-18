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
trait HasEdge extends Cloneable {
    /**
     * Get the Edge that this HasEdge has.
     */
    def edge: Edge
    
    /**
     * Do a deep clone of this HasEdge, so we can then mutate the edge without
     * affecting the original.
     *
     * We need to provide an implementation here, since the (concrete)
     * java.lang.Object implementation can't be replaced with a public abstract
     * version here. See <http://stackoverflow.com/questions/8618937/scala-
     * specifying-public-method-overriding-protected-method>
     */
    override def clone: HasEdge = throw new CloneNotSupportedException
}

/**
 * A class to say that AlleleGroups are edges.
 */
case class AlleleGroupEdge(owner: AlleleGroup) extends HasEdge {
    def edge = owner.edge
    override def clone = {
        new AlleleGroupEdge(AlleleGroup.newBuilder(owner).build)
    }
}

/**
 * A class to say that Adjacencies are edges.
 */
case class AdjacencyEdge(owner: Adjacency) extends HasEdge {
    def edge = owner.edge
    override def clone = {
        new AdjacencyEdge(Adjacency.newBuilder(owner).build)
    }
}

/**
 * A class to say that Anchors are edges.
 */
case class AnchorEdge(owner: Anchor) extends HasEdge {
    def edge = owner.edge
    override def clone = {
        new AnchorEdge(Anchor.newBuilder(owner).build)
    }
}

/**
 * A class to say that Blocks are edges.
 */
case class BlockEdge(owner: Block) extends HasEdge {
    def edge = owner.edge
    override def clone = {
        new BlockEdge(Block.newBuilder(owner).build)
    }
}
