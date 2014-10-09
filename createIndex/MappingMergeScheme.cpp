#include <stdexcept>

#include <Log.hpp>
#include <util.hpp>

#include <CreditFilter.hpp>
#include <DisambiguateFilter.hpp>

#include "MappingMergeScheme.hpp"



MappingMergeScheme::MappingMergeScheme(const FMDIndex& index, 
    const GenericBitVector& rangeVector, 
    const std::vector<std::pair<std::pair<size_t, size_t>, bool> >& rangeBases, 
    const GenericBitVector& includedPositions, size_t genome, size_t minContext, 
    size_t addContext, double multContext, double minCodingCost, bool credit, 
    std::string mapType, bool mismatch, size_t z_max, 
    std::vector<Mapping>* mappingsOut) : MergeScheme(index), threads(), 
    queue(NULL), rangeVector(rangeVector), rangeBases(rangeBases), 
    includedPositions(includedPositions), genome(genome), 
    minContext(minContext), addContext(addContext), multContext(multContext),
    minCodingCost(minCodingCost), credit(credit), mapType(mapType), 
    mismatch(mismatch), z_max(z_max), mappingsOut(mappingsOut) {
    
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
    
    // Say we need to do every contig in this genome.
    
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
        // We require different contig merge generation methods for left-right
        // schemes and for the centered schemes. Specification of mapping on
        // credit and/or allowance for inexact matching are handled within
        // these methods
        
        if(mapType == "LRexact") {
            Log::info() << "Using Left-Right exact contexts" << std::endl;
            threads.push_back(Thread(&MappingMergeScheme::generateMerges,
                this, contigsToMerge));
            
        } else {
            throw std::runtime_error(mapType + 
                " is not an implemented context type.");
          
        }
    }

    // Return a reference to the queue of merges, for our caller to do something
    // with.
    return *queue;
    
}

void MappingMergeScheme::join() {

    for(auto i = threads.begin(); 
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
        
    // Right now credit merging schemes are attempting to merge off-contig
    // positions but are otherwise behaving as expected. Throw warning and
    // don't merge rather than runtime error until this is worked out.
    
    /*if(referenceBase > index.getContigLength(referenceContig) || 
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
    }*/
    
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
    
    // Right now credit merging schemes are attempting to merge off-contig
    // positions but are otherwise behaving as expected. Throw warning and
    // don't merge rather than runtime error until this is worked out.

    
    if(queryBase > index.getContigLength(queryContig) || queryBase == 0 ||
        referenceBase > index.getContigLength(referenceContig) || 
        referenceBase == 0) {
        
        Log::error() << 
            "****WARNING: tried to merge (to) an off-contig position!****" << 
            std::endl;
    
    } else {
    
        Merge merge(queryPos, referencePos);
            
        // Send that merge to the queue.
        // Lock the queue.
        auto lock = queue->lock();
        // Spend our lock to add something to it.
        queue->enqueue(merge, lock);
        
    }
    
}

void MappingMergeScheme::generateMerges(
    ConcurrentQueue<size_t>* contigs) const {
    
    // Wait for a contig, or for there to be no more.
    auto contigLock = contigs->waitForNonemptyOrEnd();
    
    while(!contigs->isEmpty(contigLock)) {
        // We got a contig to do. Dequeue it and unlock.
        size_t queryContig = contigs->dequeue(contigLock);
        
        generateSomeMerges(queryContig);
        
        // Now wait for a new task, or for there to be no more contigs.
        contigLock = contigs->waitForNonemptyOrEnd();
        
    }
        
    // If we get here, there were no more contigs in the queue, and no
    // writers writing. Unlock it and finish the thread. TODO: Should the
    // queue just unlock if it happens to be empty?
    contigLock.unlock();
    
    // Close the output queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
}

void MappingMergeScheme::generateSomeMerges(size_t queryContig) const {
    
    
    // Set our task name to something descriptive.
    std::string taskName = "T" + std::to_string(genome) + "." + 
        std::to_string(queryContig);
        
    // How many bases have we mapped or not mapped
    size_t mappedBases = 0;
    size_t unmappedBases = 0;
    size_t creditBases = 0;
    size_t conflictedBases = 0;
    
    // Grab the contig as a string
    std::string contig = index.displayContig(queryContig);
    
    // Keep a cache of contigs by number so we don't constantly locate for
    // credit.
    std::map<size_t, std::string> contigCache;
    contigCache[queryContig] = contig;
    

    //TODO: can conflicted bases be mapped on credit?
    
    // How many positions are available to map to?
    Log::info() << taskName << " mapping " << contig.size() << 
        " bases via " << includedPositions.rank(index.getBWTLength()) <<
        " bottom-level positions" << std::endl;
    
    // We will fill these in with mappings that are to ranges, and then convert
    // the ranges to TextPositions.
    std::vector<Mapping> rightMappings;
    std::vector<Mapping> leftMappings;    
        
    if (mismatch) {
      
        Log::info() << "Using mismatch mapping with z_max " << z_max << 
            std::endl;
        
        // Map it on the right.
        rightMappings = index.misMatchMap(rangeVector, contig,
            &includedPositions, minContext, addContext, multContext, 
            minCodingCost, z_max);
        
        // Map it on the left
        leftMappings = index.misMatchMap(rangeVector, 
            reverseComplement(contig), &includedPositions, minContext, 
            addContext, multContext, minCodingCost, z_max); 
                
    } else {
        
        // We use the same function, with mismatches hardcoded to 0.
                
        // Map it on the right
        rightMappings = index.misMatchMap(rangeVector, contig,
            &includedPositions, minContext, addContext, multContext, 
            minCodingCost, 0);
        
        // Map it on the left
        leftMappings = index.misMatchMap(rangeVector, 
            reverseComplement(contig), &includedPositions, minContext, 
            addContext, multContext, minCodingCost, 0); 
    
    }
    
    // Flip the left mappings back into the original order. They should stay
    // as other-side ranges.
    std::reverse(leftMappings.begin(), leftMappings.end());
    
    for(size_t i = 0; i < leftMappings.size(); i++) {
        // Convert left and right mappings from ranges to base positions.
        
        if(leftMappings[i].isMapped()) {
            // Left is mapped, so go look up its TextPosition
            leftMappings[i].setLocation(index.getTextPosition(
                rangeBases[leftMappings[i].getRange()]));
        }

        if(rightMappings[i].isMapped()) {
            // Right is mapped, so go look up its TextPosition
            rightMappings[i].setLocation(index.getTextPosition(
                rangeBases[rightMappings[i].getRange()]));
        }
    }

    for(size_t i = 0; i < leftMappings.size(); i++) {
        // Convert all the left mapping positions to right semantics
        
        if(leftMappings[i].isMapped()) {
            // Flip anything that's mapped, using the length of the contig it
            // mapped to.
            leftMappings[i] = leftMappings[i].flip(index.getContigLength(
                leftMappings[i].getLocation().getContigNumber()));
        }
    }
    
    // We need a vector of mappings integrated from both left and right
    std::vector<Mapping> filteredMappings;
    
    if(credit) {
        // Apply a credit filter to the mappings
        filteredMappings = CreditFilter(index).apply(leftMappings,
            rightMappings);
    } else {
        // Apply a disambiguate filter to the mappings
        filteredMappings = DisambiguateFilter(index).apply(leftMappings,
            rightMappings);
    }
    
    for(size_t i = 0; i < filteredMappings.size(); i++) {
        // Now go through all the bases and make sure they have matching
        // characters before merging them. TODO: make this another filter.
        
        TextPosition candidate = filteredMappings[i].getLocation();
        
        if(!filteredMappings[i].isMapped()) {
            // Skip anything unmapped
            unmappedBases++;
            continue;
        }
                        
        // Grab the bases we are thinking of merging
        char queryBase = contig[i];
        
        // Make sure we have the contig we are mapping to. TODO: Integrate this
        // contig cache thing with the index so we can just ask the index to
        // produce the base for a TextPosition efficiently. Or just make it do
        // inverse locates fast.
        if(!contigCache.count(candidate.getContigNumber())) {
            // Grab out the other contig
            contigCache[candidate.getContigNumber()] = 
                index.displayContig(candidate.getContigNumber());
        }
        
        // Grab the base from the contig.
        char referenceBase;
        if(candidate.getStrand()) {
            // Flip around a copy and pull the complement of what's on the
            // forward strand.
            TextPosition flipped = candidate;
            flipped.flip(index.getContigLength(candidate.getContigNumber()));
            referenceBase = complement(contigCache[flipped.getContigNumber()][
                flipped.getOffset()]);
        } else {
            // Pull from the forward strand
            referenceBase = contigCache[candidate.getContigNumber()][
            candidate.getOffset()];
        }
        
        Log::debug() << queryBase << " vs. " << referenceBase << std::endl;
        
        if(queryBase == referenceBase) {
            // Only merge them if they will be the same base
            // when oriented right.
        
            // Merge coordinates are 1-based.
            generateMerge(queryContig, i + 1, 
                candidate.getContigNumber(), index.getOffset(candidate),
                candidate.getStrand());
            mappedBases++;
            
            if(mappingsOut != NULL) {
                // TODO: Assumes we have exactly one scaffold in this genome.
                
                // Offset from the start of this contig in its scaffold, and
                // save the mapping.
                (*mappingsOut)[index.getContigStart(queryContig) + i] = 
                    filteredMappings[i];
                
            }
            
        } else {
            // This base didn't map
            unmappedBases++;
        }
        
    }

    // Report that we're done.
    // Report that we're done.
    Log::info() << taskName << " finished (" << mappedBases << " mapped|" << 
        unmappedBases << " unmapped)" << std::endl;
}







