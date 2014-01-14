package edu.ucsc.soe.sequencegraph
import scala.collection.mutable.HashMap

/**
    This object handles providing sequential global IDs.
*/
object IDMaker {
    // This holds the next available ID
    var next : Long = 0
    
    /**
        Get a unique ID.
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

    SequenceGraphBuilder: a class to build up a sequence graph for a diploid
    genome.
    
    This class contains collections of all the parts needed to build a proper
    sequence graph: Loci (holding Sites and Breakpoints), Alleles, Sides, and
    SequenceGraphEdges (holding both AlleleGroups and Adjacencies).
    
    It keeps track of the last Side for each copy of each chromosome (which may
    be the same Side for both copies of a chromosome in an area with no phasing
    information). New AlleleGroups can be added to the end of each copy of each
    chromosome.
    
    Operates on a given genome/sample name, and a given reference name.

*/
class SequenceGraphBuilder(sample: String, reference: String) {
    // This holds all the Sides we have created, by ID
    val sides = HashMap.empty[Long, Side]
    
    // This holds all the SequenceGraphEdges we have created, by ID
    val edges = HashMap.empty[Long, SequenceGraphEdge]
    
    // This holds all the Alleles we have created, by ID
    val alleles = HashMap.empty[Long, Allele]
    
    // This holds all the IDs of the Sides at the ends of chromosomes, by contig
    // name and phase number (usually 0 or 1).
    val ends = HashMap.empty[(String, Int), Long]
    
    
    
    /**
        Make sure that an Allele for the given string of bases at the given Site
        Locus exists, and return its ID.
    */
    def getAllele(bases: String) : Long = {
        0
    }
    
    /**
        Attach the given AlleleGroup-carrying SequenceGraphEdge to the end of
        the given contig's given phase.
    */
    def addAlleleGroupEdge(contig: String, phase: Int, 
        alleleGroup: SequenceGraphEdge) {
        
        ends.get((contig, phase)).foreach { (end) => 
            // If we do have it already, make an edge to this
            // SequenceGraphEdge's first Side.
            val newEdge = SequenceGraphEdge.newBuilder()
                .setId(IDMaker.get())
                .setLeft(end)
                .setRight(alleleGroup.left)
                .setGenome(sample)
                // Set ploidy to 1, since we're adding to exactly 1 phase
                .setPloidy(new PloidyBounds(1, 1, null))
                .setContents(new Adjacency())
                .build()
                
            // TODO: Set up Breakpoints here, and ploidies.
            
            // Add the edge to our collection of edges
            edges(newEdge.id) = newEdge
            
        }
        
        // Put this SequenceGraphEdge's second Side as the new trailing end of
        // the chromosome.
        ends((contig, phase)) = alleleGroup.right
    }
    
    /**
        Make a new Side, with a new unique ID. Return the ID.
    */
    def makeSide() : Long = {
        // Make a new Side. TODO: Fill in position.
        val side = new Side(IDMaker.get(), new Position("somewhere", 0, Face.LEFT))
        
        // Keep it around
        sides(side.id) = side
        
        // Return its ID
        side.id
    }
    
    /**
        Make and add a new SequenceGraphEdge holding a new AlleleGroup. Takes
        the string of bases the AlleleGroup should hold, the Region of the
        reference that it corresponds to, and the ploidy of the AlleleGroup.
        
        Returns an actual SequenceGraphEdge, rather than an ID.
        
        The bases string may be null, in which case we mean to just use the
        reference sequence for the region.
    */
    def makeAlleleGroupEdge(bases: String, region: Region, ploidy: Int): 
        SequenceGraphEdge = {
        
        // Get the Allele, creating a new one if it's new.
        val allele = getAllele(bases)
        
        // Make an AlleleGroup to wrap the Allele reference
        val alleleGroup = new AlleleGroup(allele)
        
        // Make and return the actual SequenceGraphEdge
        SequenceGraphEdge.newBuilder()
            // Give it a new ID
            .setId(IDMaker.get())
            // And fresh left and right Sides
            .setLeft(makeSide())
            .setRight(makeSide())
            // Associate it with the right genome
            .setGenome(sample)
            // Set ploidy to the single integer we got
            .setPloidy(new PloidyBounds(ploidy, ploidy, null))
            // Put the AlleleGroup in it
            .setContents(alleleGroup)
            .build()
    }
}
