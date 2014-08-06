#include <stdexcept>

#include <Log.hpp>
#include <Util.h> // From libsuffixtools, for reverse_complement

#include "MappingMergeScheme.hpp"



MappingMergeScheme::MappingMergeScheme(const FMDIndex& index, 
    const BitVector& rangeVector, 
    const std::vector<std::pair<std::pair<size_t, size_t>, bool> >& rangeBases, 
    const BitVector& includedPositions, size_t genome, size_t minContext) : 
    MergeScheme(index), threads(), queue(NULL), contigsToMerge(NULL),
    rangeVector(rangeVector), rangeBases(rangeBases),
    includedPositions(includedPositions), genome(genome),
    minContext(minContext) {
    
    // Nothing to do
    
}

MappingMergeScheme::~MappingMergeScheme() {
    
    join(); // Join our threads

    if(queue != NULL) {
        // Get rid of the queue if we made one.
        delete queue;
        queue = NULL;
    }
    
    if(contigsToMerge != NULL) {
        // And the collection of contigs to merge.
        delete contigsToMerge;
        contigsToMerge = NULL;
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
    
    // Don't start more threads than we have contigs.
    size_t numThreads = std::min(MAX_THREADS, genomeContigs.second - 
        genomeContigs.first);
    
    Log::info() << "Running Mapping merge on " << numThreads << " threads" <<
        std::endl;
    
    // Make the queue of merges    
    queue = new ConcurrentQueue<Merge>(numThreads);
    
    // And one for the contigs to merge
    contigsToMerge = new ConcurrentQueue<size_t>(1);
    
    for(size_t contig = genomeContigs.first; contig < genomeContigs.second; 
        contig++) {
        
        // Put each contig in the genome into the queue of work to do.
        // TODO: Can we just use an imaginary queue of numbers in a range?
        auto lock = contigsToMerge->lock();
        contigsToMerge->enqueue(contig, lock);
    }
    
    // Say we are done writing to the queue. Good thing it doesn't have a max
    // count. Otherwise we'd need a thread to dump in all the contigs.
    auto lock = contigsToMerge->lock();
    contigsToMerge->close(lock);

    // Start up a reasonable number of threads to do the work.
    for(size_t threadID = 0; threadID < numThreads; threadID++) {
        threads.push_back(std::thread(&MappingMergeScheme::generateMerges,
                this, contigsToMerge));
    }

    // Return a reference to the queue of merges, for our caller to do something
    // with.
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
    
    if(contigsToMerge != NULL) {
        auto lock = contigsToMerge->lock();
        if(!contigsToMerge->isEmpty(lock)) {
            // Don't hold a lock and throw an exception.
            lock.unlock();
        
            // Complain our thread pool somehow broke.
            throw std::runtime_error(
                "Jobs left in queue after threads have terminated.");
        } else {
            // Everything is fine, but we still need to unlock.
            lock.unlock();
        }
    }
    
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
    
    // Where in the query do we want to come from? Always on the forward strand.
    // Correct offset to 0-based.
    TextPosition queryPos(queryContig * 2, queryBase - 1);
    
    // Where in the reference do we want to go to? May b e on the forward or
    // reverse strand. Correct text and offset for orientation, and correct
    // offset to 0-based (at the same time).
    TextPosition referencePos(referenceContig * 2 + orientation, 
        orientation ? (index.getContigLength(referenceContig) - referenceBase) :
        (referenceBase - 1));
    
    // Make a Merge between these positions.
    Merge merge(queryPos, referencePos);
        
    // Send that merge to the queue.
    // Lock the queue.
    auto lock = queue->lock();
    // Spend our lock to add something to it.
    queue->enqueue(merge, lock);
    
}

void MappingMergeScheme::generateMerges(
    ConcurrentQueue<size_t>* contigs) const {
    
    // Wait for a contig, or for there to be no more.
    auto contigLock = contigs->waitForNonemptyOrEnd();
    
    while(!contigs->isEmpty(contigLock)) {
        // We got a contig to do. Dequeue it and unlock.
        size_t queryContig = contigs->dequeue(contigLock);
    
        // Set our task name to something descriptive.
        std::string taskName = "T" + std::to_string(genome) + "." + 
            std::to_string(queryContig);
            
        // How many bases have we mapped or not mapped
        size_t mappedBases = 0;
        size_t unmappedBases = 0;
        
        // Grab the contig as a string
        std::string contig = index.displayContig(queryContig);
        
        // How many positions are available to map to?
        Log::info() << taskName << " mapping " << contig.size() << 
            " bases via " << BitVectorIterator(includedPositions).rank(
            includedPositions.getSize()) << " bottom-level positions" <<
            std::endl;
        
        // Map it on the right
        std::vector<int64_t> rightMappings = index.map(rangeVector, contig, 
            &includedPositions, minContext);
        
        // Map it on the left
        std::vector<int64_t> leftMappings = index.map(rangeVector, 
            reverseComplement(contig), &includedPositions, minContext);
        
        // Flip the left mappings back into the original order. They should stay
        // as other-side ranges.
        std::reverse(leftMappings.begin(), leftMappings.end());
          
        for(size_t i = 0; i < leftMappings.size(); i++) {
            // For each position, look at the mappings.
            
            if(leftMappings[i] != -1) {
                // We have a left mapping. Grab its base.
                auto leftBase = rangeBases[leftMappings[i]];
                
                if(rightMappings[i] != -1) {
                    // We have a right mapping. Grab its base too.
                    auto rightBase = rangeBases[rightMappings[i]];
                    
                    // Compare the position (contig, base) pairs (first) and the
                    // orientation flags (second)
                    if(leftBase.first == rightBase.first && 
                        leftBase.second != rightBase.second) {
                        
                        // These are opposite faces of the same base.
                        
                        // Produce a merge between the base we're looking at on
                        // the forward strand of this contig, and the canonical
                        // location and strand in the merged genome, accounting
                        // for orientation. TODO: Just set everything up on
                        // TextPositions or something.
                        
                        // These positions being sent are 1-based, so we have to
                        // correct i to i + 1 to get the offset of that base in
                        // the query string. Orientation is backwards to start
                        // with from our backwards right-semantics, so flip it.
                        generateMerge(queryContig, i + 1, leftBase.first.first, 
                            leftBase.first.second, !leftBase.second); 
                            
                        mappedBases++;                   
                    } else {
                        // Didn't map this one
                        unmappedBases++;
                    }
                } else {
                    // Left mapped and right didn't.
                    
                    // Do the same thing, taking the left base and merging into it.
                    // Orientation is backwards to start with from our backwards
                    // right-semantics, so flip it.
                    generateMerge(queryContig, i + 1, leftBase.first.first, 
                            leftBase.first.second, !leftBase.second);  
                            
                    mappedBases++;
                }
            } else if(rightMappings[i] != -1) {
                // Right mapped and left didn't.
                
                // We have a right mapping. Grab its base too.
                auto rightBase = rangeBases[rightMappings[i]];
                
                // Merge with the same contig and base. Leave the orientation alone
                // (since it's backwards to start with).
                generateMerge(queryContig, i + 1, rightBase.first.first, 
                            rightBase.first.second, rightBase.second); 
                
                mappedBases++; 
            } else {
                // Didn't map this one
                unmappedBases++;
            }
            
        }
        
        // Report that we're done.
        Log::info() << taskName << " finished (" << mappedBases << "|" << 
            unmappedBases << ")" << std::endl;
            
        // Now wait for a new task, or for there to be no more contigs.
        contigLock = contigs->waitForNonemptyOrEnd();
    }
    
    // If we get here, there were no more contigs in the queue, and no writers
    // writing. Unlock it and finish the thread.
    // TODO: Should the queue just unlock if it happens to be empty?
    contigLock.unlock();
    
    // Close the output queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
}










