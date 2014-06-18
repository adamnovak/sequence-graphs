#ifndef SYMMETRICMREGESCHEME_HPP
#define SYMMETRICMERGESCHEME_HPP

#include <thread>
#include <vector>

#include "MergeScheme.hpp"


/**
 * Represents the "Symmetric Merging Scheme": each genome is mapped to each
 * other genome, and bases are merged with the bases they map to.
 */
class SymmetricMergeScheme: public MergeScheme {
   
public:
    
    /**
     * Make a new SymmetricMergeScheme to merge the genomes in the given index.
     */
    SymmetricMergeScheme(const FMDIndex& index);
    
    /**
     * Get rid of a SymmetricMergeScheme (and delete its queue, if it has one).
     * If threads are running, blocks until they finish.
     */
    virtual ~SymmetricMergeScheme();
    
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
    
    /**
     * Run as a thread. Generates merges by mapping the contigs of the query
     * genome to the target genome.
     */
    virtual void generateMerges(size_t targetGenome, size_t queryGenome);
    
};

#endif
