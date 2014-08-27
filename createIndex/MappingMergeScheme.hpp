#ifndef MAPPINGMERGESCHEME_HPP
#define MAPPINGMERGESCHEME_HPP

#include "Thread.hpp"
#include <vector>

#include "MergeScheme.hpp"
#include <GenericBitVector.hpp>


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
    MappingMergeScheme(const FMDIndex& index, const GenericBitVector& rangeVector, 
        const std::vector<std::pair<std::pair<size_t, size_t>, bool> >&
        rangeBases, const GenericBitVector& includedPositions, size_t genome,
        size_t minContext = 0, size_t addContext = 0, bool credit = false, 
        std::string mapType = "LRexact",
	    bool mismatch = false, size_t z_max = 0);
    
    /**
     * Get rid of a MappingMergeScheme (and delete its queue, if it has one).
     * If threads are running, blocks until they finish.
     */
    virtual ~MappingMergeScheme();
    
    /**
     * Create and return a queue of merges, and start feeding merges into it
     * from some other thread(s). The writers on the queue will be known to it,
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

    // How many worker threads should be started, maximum, to produce merges
    // from contigs? TODO: Magically know a good answer for all systems.
    static const size_t MAX_THREADS = 32;

    // Holds all the threads that are generating merges.
    std::vector<Thread> threads;
    
    // Holds a pointer to a ConcurrentQueue, so we can create one and then
    // destroy it only when we get destroyed.
    ConcurrentQueue<Merge>* queue;
    
    // Holds a ConcurrentQueue of all the contig IDs that need to be processed.
    // We fill this up with contig IDs, and our merge threads read from it, so
    // we can pool them instead of spawning about a thousand of them. It will
    // have one writer, which is done pretty much as soon as the queue starts
    // getting used.
    ConcurrentQueue<size_t>* contigsToMerge;
    
    // Holds the bit vector marking out the BWT ranges that belong to higher-
    // level positions.
    const GenericBitVector& rangeVector;
    
    // Holds the canonicalized position to merge into for a forward mapping to
    // each range. Positions are 1-based.
    const std::vector<std::pair<std::pair<size_t, size_t>, bool> >& rangeBases; 
        
    // Holds a mask; we should only consider positions with 1s here as actually
    // being in the bottom level.
    const GenericBitVector& includedPositions;
    
    // Holds the number of the genome we are going to map.
    size_t genome;
    
    // Minimum amount of context that is allowed to motivate a merge.
    size_t minContext;
    
    // Minimum distance you have to go out past where a context is unique in
    // order to actually map.
    size_t addContext;
    
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
     * genome; left-right contexts
     */
    virtual void generateMerges(ConcurrentQueue<size_t>* contigs) const;
    
    /**
     * Generate laft-right merges from one particular contig.
     */
    virtual void generateSomeMerges(size_t queryContig) const;

    /**
     * Run as a thread. Generates merges by mapping a query contig to the target
     * genome; centered contexts
     */
    virtual void CgenerateMerges(ConcurrentQueue<size_t>* contigs) const;
    
    /**
     * Generate centered-context merges for one particular contig.
     */
    virtual void CgenerateSomeMerges(size_t queryContig) const;
    
};

#endif
