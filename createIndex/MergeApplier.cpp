#include <iostream>

#include <Log.hpp>

#include "MergeApplier.hpp"

MergeApplier::MergeApplier(const FMDIndex& index,
    ConcurrentQueue<Merge>& source, stPinchThreadSet* target): index(index), 
    source(source), target(target), thread(&MergeApplier::run, this) {
    
    // Already started running. See <http://stackoverflow.com/a/10673671/402891>
    
}

void MergeApplier::join() {
    thread.join();
}

void MergeApplier::run() {
    // OK, do the actual merging.
    
    while(true) {
        // Lock the queue when either there's something in it or all its writers
        // have closed it. If neither happens we just wait here forever.
        auto lock = source.waitForNonemptyOrEnd();
        
        if(source.isEmpty(lock)) {
            // We've finished our job.
            return;
        }
        
        // If we get here, there's actual work to do. Trade our lock for a value
        // to work on.
        Merge merge = source.dequeue(lock);
        
        // Now actually apply the merge.
        
        // Unpack the first TextPosition to be merged
        size_t firstContigNumber = index.getContigNumber(merge.first);
        size_t firstStrand = index.getStrand(merge.first);
        size_t firstOffset = index.getOffset(merge.first);
        
        // And the second
        size_t secondContigNumber = index.getContigNumber(merge.second);
        size_t secondStrand = index.getStrand(merge.second);
        size_t secondOffset = index.getOffset(merge.second);
        
        // What orientation should we use for the second strand, given that we
        // are pinching against the first strand in orientation 1 (reverse)?
        bool orientation = firstStrand == secondStrand;
        
        // Grab the first pinch thread
        stPinchThread* firstThread = stPinchThreadSet_getThread(target,
            firstContigNumber);
            
        // And the second
        stPinchThread* secondThread = stPinchThreadSet_getThread(target,
            secondContigNumber);
            
        // Log the pinch if applicable.
        Log::trace() << "\tPinching #" << firstContigNumber << ":" <<
            firstOffset << " strand " << firstStrand << " and #" << 
            secondContigNumber << ":" << secondOffset << " strand " << 
            secondStrand << " (orientation: " << orientation << ")" <<
            std::endl;
        
        // Perform the pinch
        stPinchThread_pinch(firstThread, secondThread, firstOffset,
            secondOffset, 1, orientation);
        
        
        
        
        
    }
    
}


