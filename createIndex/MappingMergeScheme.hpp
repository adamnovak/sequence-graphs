#ifndef MAPPINGMERGESCHEME_HPP
#define MAPPINGMERGESCHEME_HPP

#include <thread>
#include <vector>

#include "MergeScheme.hpp"


/**
 * Represents the "Mapping Merging Scheme", a component of the greedy merging
 * scheme: every contig in a genome is mapped to the given merged reference
 * structure.
 */
class MappingMergeScheme: public MergeScheme {
   
public:
    
    /**
     * Make a new MappingMergeScheme, which takes the index defined by the given
     * range vector and vector of per-position 1-based canonicalized bases to
     * merge with, and maps each contig in the given genome to it. Only maps on
     * contexts provided by the positions marked in includedPositions, and
     * requires at least the given minimum number of bases of context to map.
     *
     * TODO: Make this function take a less absurd number of things.
     */
    MappingMergeScheme(const FMDIndex& index, const BitVector& rangeVector, 
        const std::vector<std::pair<std::pair<size_t, size_t>, bool> >&
        rangeBases, const BitVector& includedPositions, size_t genome,
        size_t minContext = 0, bool credit = false, std::string mapType = "LRexact",
	bool mismatch = false, size_t z_max = 0);
    
    /**
     * Get rid of a MappingMergeScheme (and delete its queue, if it has one).
     * If threads are running, blocks until they finish.
     */
    virtual ~MappingMergeScheme();
    
    /**
     * Create and return a queue of merges, and start feeding merges into it
     * from some other thread(s). The writers on the queue must be know to it,
     * so that the queue will know when all merges have been written.
     *
     * May only be called once.
     * 
     * Returns a reference to the queue, which will live as long as this object
     * does.
     */
    virtual ConcurrentQueue<Merge>& run() override;
    
    /**
     * Wait for all the merge-producing threads to finish. Obviously you
     * shouldn't call this unless you've finished reading the queue from run and
     * know no more merges will be generated.
     */
    virtual void join() override;
    
protected:

    // Holds all the threads that are generating merges.
    std::vector<std::thread> threads;
    
    // Holds a pointer to a ConcurrentQueue, so we can create one and then
    // destroy it only when we get destroyed.
    ConcurrentQueue<Merge>* queue;
    
    // Holds the bit vector marking out the BWT ranges that belong to higher-
    // level positions.
    const BitVector& rangeVector;
    
    // Holds the canonicalized position to merge into for a forward mapping to
    // each range. Positions are 1-based.
    const std::vector<std::pair<std::pair<size_t, size_t>, bool> >& rangeBases; 
        
    // Holds a mask; we should only consider positions with 1s here as actually
    // being in the bottom level.
    const BitVector& includedPositions;
    
    // Holds the number of the genome we are going to map.
    size_t genome;
    
    // Minimum amount of context that is allowed to motivate a merge.
    size_t minContext;
    
    // Flag whether to use mapping on credit scheme   
    bool credit;
    
    // Type of context used
    std::string mapType;
    
    bool mismatch;
    
    size_t z_max;
    
    /**
     * Create a Merge between two positions and enqueue it. Positions are
     * 1-based.
     */
    void generateMerge(size_t queryContig, size_t queryBase, 
        size_t referenceContig, size_t referenceBase, bool orientation) const;
    
    /**
     * Run as a thread. Generates merges by mapping a query contig to the target
     * genome; left-right exact contexts
     */
    virtual void generateMerges(size_t queryContig) const;
    
    /**
     * Run as a thread. Generates merges by mapping a query contig to the target
     * genome; centered contexts
     */
    virtual void CgenerateMerges(size_t queryContig) const;
    
};

#endif
