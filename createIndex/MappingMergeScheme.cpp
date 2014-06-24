#include <stdexcept>

#include <Log.hpp>
#include <Util.h> // From libsuffixtools, for reverse_complement

#include "MappingMergeScheme.hpp"



MappingMergeScheme::MappingMergeScheme(const FMDIndex& index, 
    const BitVector& rangeVector, 
    const std::vector<std::pair<std::pair<size_t, size_t>, bool> >& rangeBases, 
    const BitVector& includedPositions, size_t genome, size_t minContext) : 
    MergeScheme(index), threads(), queue(NULL), rangeVector(rangeVector), 
    rangeBases(rangeBases), includedPositions(includedPositions), 
    genome(genome), minContext(minContext) {
    
    // Nothing to do
    
}

MappingMergeScheme::~MappingMergeScheme() {
    
    join(); // Join our threads

    if(queue != NULL) {
        // Get rid of the queue if we made one.
        delete queue;
        queue = NULL;
    }
    
}

ConcurrentQueue<Merge>& MappingMergeScheme::run() {
    
    if(queue != NULL) {
        // Don't let people call this twice.
        throw std::runtime_error(
            "Called run() twice on a MappingMergeScheme.");
    }
    
    // Grab the limits of the contig range belonging to this genome.
    auto genomeContigs = index.getGenomeContigs(genome);
    
    // Figure out how many threads to run. Let's do one per contig in this
    // genome.
    size_t threadCount = genomeContigs.second - genomeContigs.first;
    
    Log::info() << "Running Mapping merge on " << threadCount << " threads" <<
        std::endl;
    
    // Make the queue    
    queue = new ConcurrentQueue<Merge>(threadCount);
    
    for(size_t contig = genomeContigs.first; contig < genomeContigs.second; 
        contig++) {
        
        // Start up a thread for every contig in this genome.
        threads.push_back(std::thread(&MappingMergeScheme::generateMerges,
                this, contig));
    }

    // Return a reference to the queue.
    return *queue;
    
}

void MappingMergeScheme::join() {

    for(std::vector<std::thread>::iterator i = threads.begin(); 
        i != threads.end(); ++i) {
    
        // Join each thread.
        (*i).join();    
        
    }
    
    // Probably not a good idea to join the same threads twice.
    threads.clear();
    
}

void MappingMergeScheme::generateMerge(size_t queryContig, size_t queryBase, 
    size_t referenceContig, size_t referenceBase, bool orientation) const {
    
    if(referenceBase > index.getContigLength(referenceContig) || 
        referenceBase == 0) {
        
        // Complain that we're trying to talk about an off-the-contig position.
        throw std::runtime_error("Reference base " + 
            std::to_string(referenceBase) + 
            " on contig " + std::to_string(referenceContig) + 
            " is beyond 1-based bounds of length " + 
            std::to_string(index.getContigLength(referenceContig)));
    }
    
    if(queryBase > index.getContigLength(queryContig) || 
        queryBase == 0) {
        
        // Complain that we're trying to talk about an off-the-contig position.
        throw std::runtime_error("Query base " + std::to_string(queryBase) + 
            " on contig " + std::to_string(queryContig) + 
            " is beyond 1-based bounds of length " + 
            std::to_string(index.getContigLength(queryContig)));
    }
    
    // Make a Merge between the specified query base and the specified reference
    // base. Account for orientation and change to 0-based coordinates.
    Merge merge(TextPosition(queryContig * 2, queryBase - 1), 
        TextPosition(referenceContig * 2 + orientation, referenceBase - 1));
        
    // Send that merge to the queue.
    // Lock the queue.
    auto lock = queue->lock();
    // Spend our lock to add something to it.
    queue->enqueue(merge, lock);
    
}

void MappingMergeScheme::generateMerges(size_t queryContig) const {
    
    // Spot check all of the ranges
    for(auto i : rangeBases) {
        // Grab the contig number
        size_t rangeContig = i.first.first;
        // And the 1-based offset into it
        size_t rangeOffset = i.first.second;
        
        if(rangeOffset == 0 || 
            rangeOffset > index.getContigLength(rangeContig)) {
            
            // We found a bad one. Complain.
            throw std::runtime_error("Out of bounds canonical base!");
        }
    }
    
    // What's our thread name?
    std::string threadName = "T" + std::to_string(queryContig) + "->" + 
        std::to_string(genome);
    
    // Grab the contig as a string
    std::string contig = index.displayContig(queryContig);
    
    // Map it on the right
    std::vector<int64_t> rightMappings = index.map(rangeVector, contig, 
        &includedPositions, minContext);
    
    // Map it on the left
    std::vector<int64_t> leftMappings = index.map(rangeVector, 
        reverseComplement(contig), &includedPositions, minContext);
        
    // Flip the left mappings back into the original order. They should stay as
    // other-side ranges.
    std::reverse(leftMappings.begin(), leftMappings.end());
    
    for(size_t i = 0; i < leftMappings.size(); i++) {
        // For each position, look at the mappings.
        
        if(leftMappings[i] != -1) {
            // We have a left mapping. Grab its base.
            auto leftBase = rangeBases[leftMappings[i]];
            
            if(rightMappings[i] != -1) {
                // We have a right mapping. Grab its base too.
                auto rightBase = rangeBases[rightMappings[i]];
                
                if(leftBase.first == rightBase.first && 
                    leftBase.second != rightBase.second) {
                    
                    // These are opposite faces of the same base.
                    
                    // Produce a merge between the base we're looking at on the
                    // forward strand of this contig, and the canonical location
                    // and strand in the merged genome, accounting for
                    // orientation. TODO: Just set everything up on
                    // TextPositions or something.
                    
                    // These positions being sent are 1-based, so we have to
                    // correct i to i + 1 to get the offset of that base in the
                    // query string. Orientation is backwards to start with from
                    // our backwards right-semantics, so flip it.
                    generateMerge(queryContig, i + 1, leftBase.first.first, 
                        leftBase.first.second, !leftBase.second);                    
                }
            } else {
                // Left mapped and right didn't.
                
                // Do the same thing, taking the left base and merging into it.
                // Orientation is backwards to start with from our backwards
                // right-semantics, so flip it.
                generateMerge(queryContig, i + 1, leftBase.first.first, 
                        leftBase.first.second, !leftBase.second);  
            }
        } else if(rightMappings[i] != -1) {
            // Right mapped and left didn't.
            
            // We have a right mapping. Grab its base too.
            auto rightBase = rangeBases[rightMappings[i]];
            
            // Merge with the same contig and base. Leave the orientation alone
            // (since it's backwards to start with).
            generateMerge(queryContig, i + 1, rightBase.first.first, 
                        rightBase.first.second, rightBase.second);  
        }
        
    }
    
    // Close the queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
    // Report that we're done.
    Log::info() << threadName << " finished" << std::endl;
}










