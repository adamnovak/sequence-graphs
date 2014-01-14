package edu.ucsc.genome
import scala.collection.mutable.HashMap
import scala.collection.mutable.MutableList

/**
 * This object handles providing sequential global IDs.
 */
object IDMaker {
    // This holds the next available ID
    var next : Long = 0
    
    /**
     * Get a unique ID.
     */
    def get() : Long = {
        // Make an ID to return
        val id = next
        // Advance so we generate a different ID next time.
        next += 1
        // Return the ID we generated
        id
    }
}

/**
 *
 * SequenceGraphBuilder: a class to build up a sequence graph for a diploid
 * genome.
 * 
 * This class contains collections of all the parts needed to build a proper
 * sequence graph.
 * 
 * It keeps track of the last Side for each copy of each chromosome (which may
 * be the same Side for both copies of a chromosome in an area with no phasing
 * information). New AlleleGroups can be added to the end of each copy of each
 * chromosome.
 * 
 * Operates on a given genome/sample name, and a given reference name.
 *
 */
class SequenceGraphBuilder(sample: String, reference: String) {
    // This holds all the Sides we have created, by ID
    val sides = HashMap.empty[Long, Side]
    
    // This holds all the AlleleGroups we have created.
    val alleleGroups = MutableList[AlleleGroup]()
    
    // This holds all the Adjacencies we have created.
    val adjacencies = MutableList[Adjacency]()
    
    // This holds all the Anchors we have created.
    val anchors = MutableList[Anchor]()
    
    // This holds all the IDs of the Sides at the ends of chromosomes, by contig
    // name and phase number (usually 0 or 1).
    val ends = HashMap.empty[(String, Int), Long]
    
    /**
     * Attach the given AlleleGroup to the end of the given contig's given
     * phase.
     */
    def addAlleleGroup(contig: String, phase: Int, alleleGroup: AlleleGroup) = {
        
        ends.get((contig, phase)).foreach { (end) => 
            // If we do have something at the end of this phase of this contig
            // already, make an Adjacency to this AlleleGroup's first Side.
            val newAdjacency = Adjacency.newBuilder()
                // Attach the Edge
                .setEdge(new Edge(IDMaker.get(), end, alleleGroup.edge.left))
                // Set ploidy to exactly 1, which can be a null PloidyBounds.
                .setPloidy(null)
                // Attach to our genome
                .setGenome(sample)
                .build()
                
            // Add the Adjacency to our collection
            adjacencies += newAdjacency
            
        }
        
        // Put this AlleleGroup's second Side as the new trailing end of
        // the chromosome.
        ends((contig, phase)) = alleleGroup.edge.right
        
        // Remember the AlleleGroup in our list of AlleleGroups
        alleleGroups += alleleGroup
    }
    
    /**
     * Remember the given Side by its ID. Can be called repeatedly on the same
     * Side.
     */
    def addSide(side: Side) = {
        // Put the Side in the HashMap
        sides(side.id) = side
    }
    
    /**
     * Create and return (but do not remember) a Side corresponding to the 5'
     * face of the next unaccounted-for base of the given phase of the given
     * contig.
     *
     * If there is nothing at the end of that phase of that contig, or if the
     * thing present isn't a right side of a base that's actually on that
     * contig, produces Nothing.
     */
    def getNextSide(contig: String, phase: Int) : Option[Side] = {
        // Go get the side ID of the end of the contig, and then the Side
        // itself, and then work on that.
        ends.get((contig, phase)).flatMap(sides.get).flatMap({ (end : Side) =>
            if(end.position.face == Face.RIGHT &&
                end.position.contig == contig) {
                
                // The Side already there looks like it belongs there; we can
                // just go 1 base right and present the left Side.
                Some(new Side(IDMaker.get(), new Position(contig, 
                    end.position.base + 1, Face.LEFT)))
                
            } else {
                // The Side there doesn't look reasonable, we can't just extend
                // it right arbitrarily.
                None
            }
        })
    }
}








