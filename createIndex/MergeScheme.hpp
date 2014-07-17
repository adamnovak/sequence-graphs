#ifndef MERGESCHEME_HPP
#define MERGESCHEME_HPP

#include <FMDIndex.hpp>

#include "ConcurrentQueue.hpp"
#include "Merge.hpp"

/**
 * Represents a merging scheme which starts a bunch of threads and dumps Merges
 * into a ConcurrentQueue. Not all merging schemes will fit this base class;
 * some need to control the outer loop and build multiple indexes one after the
 * other.
 */
class MergeScheme {

public:

    /**
     * Make a new MergeScheme that will produce merges on the given index.
     */
    MergeScheme(const FMDIndex& index);

    /** 
     * We have a default virtual destructor, which we are going to need if we're
     * going to ever use polymorphism.
     */
    virtual ~MergeScheme() = default;
    
    /**
     * MergeSchemes can be move-constructed.
     */
    MergeScheme(MergeScheme&& other) = default;
    
    /**
     * MergeSchemes can be move-assigned.
     */
    MergeScheme& operator=(MergeScheme&& other) = default;
    
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
    virtual ConcurrentQueue<Merge>& run() = 0;
    
    /**
     * Wait for all the merge-producing threads to finish. Obviously you
     * shouldn't call this unless you've finished reading the queue from run and
     * know no more merges will be generated.
     */
    virtual void join() = 0;
    
protected:

    // Holds the FMDIndex we're using to look at the low-level sequences.
    const FMDIndex& index;
    
private:
    /**
     * MergeSchemes cannot be copy-constructed.
     */
    MergeScheme(const MergeScheme& other) = delete;
    
    /**
     * MergeSchemes cannot be copy-assigned.
     */
    MergeScheme& operator=(const MergeScheme& other) = delete;
    
            

};

#endif
