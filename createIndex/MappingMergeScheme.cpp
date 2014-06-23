#include "MappingMergeScheme.hpp"

#include <Log.hpp>
#include <Util.h> // From libsuffixtools, for reverse_complement

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
    
    // Make a Merge between the specified query base and the specified reference
    // base. Account for orientation.
    Merge merge(TextPosition(queryContig * 2, queryBase), 
        TextPosition(referenceContig * 2 + orientation, referenceBase));
        
    // Send that merge to the queue.
    // Lock the queue.
    auto lock = queue->lock();
    // Spend our lock to add something to it.
    queue->enqueue(merge, lock);
    
}

void MappingMergeScheme::generateMerges(size_t queryContig) const {
    
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
                    generateMerge(queryContig, i, leftBase.first.first, 
                        leftBase.first.second, leftBase.second);                    
                }
            } else {
                // Left mapped and right didn't.
                
                // Do the same thing, taking the left base and merging into it.
                generateMerge(queryContig, i, leftBase.first.first, 
                        leftBase.first.second, leftBase.second);  
            }
        } else if(rightMappings[i] != -1) {
            // Right mapped and left didn't.
            
            // We have a right mapping. Grab its base too.
            auto rightBase = rangeBases[rightMappings[i]];
            
            // Merge with the same contig and base, but flip the orientation,
            // since we want left=false=default semantics.
            generateMerge(queryContig, i, rightBase.first.first, 
                        rightBase.first.second, !rightBase.second);  
        }
        
    }
    
    // Close the queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
    // Report that we're done.
    Log::info() << threadName << " finished" << std::endl;
}










