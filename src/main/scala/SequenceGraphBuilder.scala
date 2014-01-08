package edu.ucsc.soe.sequencegraph
import scala.collection.mutable.HashMap

/**
    This objecthandle sproviding sequential global IDs.
*/
object IDMaker {
    // This holds the next available ID
    var next = 0
    
    /**
        Get a unique string ID.
    */
    def get() : String = {
        // Make an ID string to return
        val id = next.toString()
        // Advance so we ghenerate a different ID next time.
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

*/
class SequenceGraphBuilder(sample: String) {
    // This holds all the Sites and Breakpoints we have created inside Loci, by
    // ID
    val loci = HashMap.empty[String, Locus]
    
    // This holds all the Sides we have created, by ID
    val sides = HashMap.empty[String, Side]
    
    // This holds all the Alleles we have created, by ID
    val alleles = HashMap.empty[String, Allele]
    
    // This holds all the IDs of the Sides at the ends of chromosomes, by contig
    // name and phase number (usually 0 or 1).
    val ends = HashMap.empty[(String, Int), String]
    
    /**
        Make sure that a Site for the given region in the given reference
        exists, and return the ID of its Locus.
    */
    def getSite(region: Region, reference: String): String = {
        // We cheat a bit by making Site Locus IDs just be the stringified
        // contents of the Site.
        val id = "%s:%s:%d-%d".format(reference, region.contig, region.start,
            region.end)
        
        if(!loci.contains(id)) {
            // We need to add this Site first
            loci(id) = new Locus(id, new Site(region, reference))
        }
        
        // Return the ID of the Site Locus which is now guaranteed to exist.
        id
    }
    
    /**
        Make sure that an Allele for the given string of bases at the given Site
        Locus exists, and return its ID.
    */
    def getAllele(locus: String, bases: String) : String = {
        // We cheat a bit by making Allele IDs just be the stringified
        // contents of the Allele.
        val id = "%s=%s".format(locus, bases)
        
        if(!alleles.contains(id)) {
            // We need to add this Allele first
            alleles(id) = new Allele(id, locus, bases)
        }
        
        // Return the ID of the Allele which is now guaranteed to exist.
        id
    }
    
    /**
        Attach the given AlleleGroup-carrying SequenceGraphEdge to the end of
        the given contig's given phase.
    */
    def addAlleleGroup(contig: String, phase: Int, 
        alleleGroup: SequenceGraphEdge) {
        
        ends.get((contig, phase)).foreach { (end) => 
            // If we do have it already, add an edge to this SequenceGraphEdge's
            // first Side.
            val newEdge = SequenceGraphEdge.newBuilder()
                .setId(IDMaker.get())
                .setLeft(end)
                .setRight(IDMaker.get())
                // Set ploidy to 1, since we're adding to exactly 1 phase
                .setPloidy(new PloidyBounds(1, 1))
                .setContents(new Adjacency())
            // TODO: Set up Breakpoints here, and ploidies.
        }
        
        // Put this SequenceGraphEdge's second Side as the new trailing end of
        // the chromosome.
        ends((contig, phase)) = alleleGroup.right
    }
}
