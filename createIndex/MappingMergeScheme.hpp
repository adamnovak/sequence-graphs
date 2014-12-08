#ifndef MAPPINGMERGESCHEME_HPP
#define MAPPINGMERGESCHEME_HPP

#include "Thread.hpp"
#include <vector>

#include "MergeScheme.hpp"
#include <GenericBitVector.hpp>
#include <MappingScheme.hpp>


/**
 * Represents the "Mapping Merging Scheme", a component of the greedy merging
 * scheme: every contig in a genome is mapped to the given merged reference
 * structure.
 */
class MappingMergeScheme: public MergeScheme {
   
public:
    
    /**
     * Make a new MappingMergeScheme, which maps the given genome from the given
     * index using the given mapping scheme.
     */
    MappingMergeScheme(const FMDIndex& index,
        const MappingScheme* mappingScheme, size_t genome);
    
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
    static const size_t MAX_THREADS;

    // Keep the index so we can pull the contigs from it.
    const FMDIndex& index;

    // Holds the number of the genome we are going to map.
    size_t genome;
    
    // Holds a ConcurrentQueue of all the contig IDs that need to be processed.
    // We fill this up with contig IDs, and our merge threads read from it, so
    // we can pool them instead of spawning about a thousand of them. It will
    // have one writer, which is done pretty much as soon as the queue starts
    // getting used.
    ConcurrentQueue<size_t>* contigsToMerge;

    // Holds all the threads that are generating merges.
    std::vector<Thread> threads;
    
    // Holds a pointer to a ConcurrentQueue, so we can create one and then
    // destroy it only when we get destroyed.
    ConcurrentQueue<Merge>* queue;
    
    // Holds the mapping scheme we will use to map.
    const MappingScheme* mappingScheme;

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
     * Generate left-right merges from one particular contig.
     */
    virtual void generateSomeMerges(size_t queryContig) const;

};

#endif
