package edu.ucsc.soe.sequencegraph
import scala.collection.mutable.HashMap

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
}
