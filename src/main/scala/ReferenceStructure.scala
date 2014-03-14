package edu.ucsc.genome

/**
 * Represents a Reference Structure: a phased or unphased sequence graph, with a
 * rule for mapping to it. Does not necessarily keep the SequenceGraph around.
 *
 * Internally, has two implementations: one for the bottom-level collection of
 * haplotypes, and one for a higher-level graph that merges nodes from the level
 * below.
 */
trait ReferenceStructure {
    /**
     * Map a string on both sides. Returns a sequence of Positions to which each
     * base maps, with None for bases that don't map. All the Position UUIDs
     * returned, which are (contig, base) pairs, will be specific to this level.
     */
    def map(pattern: String): Seq[Option[Position]]
    
}

/**
 * A ReferenceStructure for phased sequence graphs (i.e. string haplotypes).
 * Backed by an FMDIndex.
 *
 * Takes a string basename for the index to load, as created with RLCSABuilder.
 * Positions are all on contigs in that index.
 */
class StringReferenceStructure(basename: String) extends ReferenceStructure {
    
    // Make our FMDIndex
    val index = new FMDIndex(basename)

    // Map with our index
    def map(pattern: String): Seq[Option[Position]] = index.map(pattern)

}

/**
 * A ReferenceStructure which has some things in a lower-level
 * ReferenceStructure collapsed together.
 */
abstract class CollapsedReferenceStructure(base: ReferenceStructure) 
    extends ReferenceStructure {

    // We have an IntervalMap of intervals we merge to what we merge them into.
    val intervals = new IntervalMap[Position]
}
