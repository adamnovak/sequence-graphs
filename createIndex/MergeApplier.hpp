#ifndef MERGEAPPLIER_HPP
#define MERGEAPPLIER_HPP

#include <stPinchGraphs.h>

#include <FMDIndex.hpp>

#include "ConcurrentQueue.hpp"
#include "Merge.hpp"
#include "Thread.hpp"

/**
 * A class which reads in from a ConcurrentQueue of Merges and applies them all
 * to an stPinchGraph.
 */
class MergeApplier {

public:
    /**
     * Make a new MergeApplier to, using the given index, apply merges from the
     * given queue to the given pinch graph. Automatically starts running. The
     * ConcurrentQueue must have been initialized with some number of writers.
     */
    MergeApplier(const FMDIndex& index, ConcurrentQueue<Merge>& source,
        stPinchThreadSet* target);
    
    /**
     * Wait for the merge applier to finish its work.
     */
    void join();
    
protected:
    // Keep a reference to the index we'll use to turn TextPositions into
    // coordinates on the pinch graph.
    const FMDIndex& index;
    
    // Keep around the queue where merges come from.
    ConcurrentQueue<Merge>& source;
    
    // Keep around a pointer to the graph to apply the merges to.
    stPinchThreadSet* target;
    
    // Keep around a thread that runs to do the actual applying.
    Thread thread;
    
    /**
     * Thread that does the actual work. Gets to own the target graph while
     * running.
     */
    void run();
    
};

#endif
