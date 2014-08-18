#ifndef OVERLAPMERGESCHEME_HPP
#define OVERLAPMERGESCHEME_HPP

#include <thread>
#include <vector>

#include "MergeScheme.hpp"


/**
 * Represents the "Overlap Merging Scheme": each genome is mapped to each
 * other genome, and bases are merged with the bases they map to.
 */
class OverlapMergeScheme: public MergeScheme {
   
public:
    
    /**
     * Make a new OverlapMergeScheme to merge the genomes in the given index. If
     * minContext is specified, ignores merges motivated by fewer than that
     * number of bases of context, even if there is an unambiguous mapping.
     */
    OverlapMergeScheme(const FMDIndex& index, size_t minContext = 0);
    
    /**
     * Get rid of a OverlapMergeScheme (and delete its queue, if it has one).
     * If threads are running, blocks until they finish.
     */
    virtual ~OverlapMergeScheme();
    
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
    
    // Minimum amount of context that is allowed to motivate a merge.
    size_t minContext;
    
    /**
     * Run as a thread. Generates merges by mapping the contigs of the query
     * genome to the target genome.
     */
    virtual void generateMerges(size_t targetGenome, size_t queryGenome);
    
};

#endif
