#include <stdexcept>

#include <Log.hpp>
#include <util.hpp>

#include <CreditFilter2.hpp>
#include <DisambiguateFilter.hpp>

#include "MappingMergeScheme.hpp"

const size_t MappingMergeScheme::MAX_THREADS = 32;

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

    if(mapType == "LRexact") {
        // Give a message saying we're doing LR matching (which, contrary to its
        // "exact" name, supports mismatches)
        Log::info() << "Using Left-Right exact contexts" << std::endl;
    } else if(mapType == "natural") {
        // Give a message saying we're doing Benedict's new natural mapping.
        Log::info() << "Using natural contexts" << std::endl;
        if(mismatch > 0) {
            // Didn't write this code.
            throw std::runtime_error("Natural mapping with mismatches is "
                "probably not anywhere near efficient and thus not "
                "implemented");
        }
    } else {
        throw std::runtime_error(mapType + 
            " is not an implemented context type.");
    }

    // Start up a reasonable number of threads to do the work.
    for(size_t threadID = 0; threadID < numThreads; threadID++) {
        // All the mapTypes can use the same method here, since there's lots of
        // work scheduling and setup logic in common.
        
        // Start up a thread.
        threads.push_back(Thread(&MappingMergeScheme::generateMerges,
            this, contigsToMerge));
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

/**
 * Get the character at the given TextPosition, using the contig cache.
 *
 * TODO: Move into FMDIndex now that it has that cache and can support this 
 * somewhat efficiently.
 */
char getCharacter(const FMDIndex& index, TextPosition position) {
    
    if(position.getStrand()) {
        // On the reverse strand.
        
        // Flip our copy of the position.
        position.flip(index.getContigLength(position.getContigNumber()));
        
        // Get the character on the other strand and complement it.
        return complement(getCharacter(index, position));
    } else {
        // On the forward strand
        
        // Grab the reference contig in the cache
        const std::string& referenceContig = 
            index.displayContigCached(position.getContigNumber());
            
        // Pull out the correct character.
        return referenceContig[position.getOffset()];
    }
}

void MappingMergeScheme::generateSomeMerges(size_t queryContig) const {
    
    
    // Set our task name to something descriptive.
    std::string taskName = "T" + std::to_string(genome) + "." + 
        std::to_string(queryContig);
        
    // Grab the contig as a string
    std::string contig = index.displayContig(queryContig);
    
    // How many positions are available to map to?
    Log::info() << taskName << " mapping " << contig.size() << 
        " bases via " << includedPositions.rank(index.getBWTLength()) <<
        " bottom-level positions" << std::endl;
    
    // How many bases have we mapped or not mapped (credit or not).
    size_t mappedBases = 0;
    size_t unmappedBases = 0;
    
    // How many bases are mapped on credit?
    size_t creditBases = 0;
    // How many bases have conflicted credit?
    size_t conflictedCredit = 0;
    
    if(mapType == "LRexact") {
    
        
        // How many are mapped on both the left and the right?
        size_t leftRightMappedBases = 0;
        
        // We will fill these in with mappings that are to ranges, and then
        // convert the ranges to TextPositions.
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
                // Flip anything that's mapped, using the length of the contig
                // it mapped to.
                leftMappings[i] = leftMappings[i].flip(index.getContigLength(
                    leftMappings[i].getLocation().getContigNumber()));
            }
        }
        
        // We need a vector of mappings integrated from both left and right
        std::vector<Mapping> filteredMappings;
        
        if(credit) {
            // Apply a credit filter to the mappings
            filteredMappings = CreditFilter2(index, rangeVector, 
                mismatch? z_max : 0,  &includedPositions).apply(leftMappings, 
                rightMappings, index.displayContigCached(queryContig));
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
            
            // Make sure we have the contig we are mapping to.
            const std::string& referenceContig = index.displayContigCached(
                candidate.getContigNumber());
            
            // Grab the base from the contig.
            char referenceBase;
            if(candidate.getStrand()) {
                // Flip around a copy of the position and pull the complement of
                // what's on the forward strand, so we don't go and complement
                // the whole contig.
                TextPosition flipped = candidate;
                flipped.flip(index.getContigLength(
                    candidate.getContigNumber()));
                
                referenceBase = complement(
                    referenceContig[flipped.getOffset()]);
            } else {
                // Pull from the forward strand
                referenceBase = referenceContig[candidate.getOffset()];
            }
            
            Log::trace() << queryBase << " vs. " << referenceBase << std::endl;
            
            if(queryBase == referenceBase) {
                // Only merge them if they will be the same base
                // when oriented right.
            
                // Merge coordinates are 1-based.
                generateMerge(queryContig, i + 1, 
                    candidate.getContigNumber(), index.getOffset(candidate),
                    candidate.getStrand());
                mappedBases++;
                
                if(leftMappings[i].isMapped() && rightMappings[i].isMapped()) {
                    // This base was left-mapped and right-mapped (and also not
                    // mapped on credit).
                    leftRightMappedBases++;
                }
                
                if(mappingsOut != NULL) {
                    // TODO: Assumes we have exactly one scaffold in this
                    // genome.
                    
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
        Log::info() << taskName << " finished (" << mappedBases << " mapped|" << 
            unmappedBases << " unmapped|" << leftRightMappedBases <<
            " LR-mapped)" << std::endl;
    } else if(mapType == "natural") {
        // Map using the natural context scheme: get matchings from all the
        // unique-in-the-reference strings that overlap you.
        
        // Call into the index. TODO: pass a parameters struct of some type.
        // TODO: It would make sense to use a DisambiguateFilter or some such,
        // but with this algorithm that's sort of integrated into the mapping.
        // TODO: This doesn't even use ranges.
        std::vector<Mapping> naturalMappings = index.naturalMap(contig,
            &includedPositions, minContext);
        
        for(size_t i = 0; i < naturalMappings.size(); i++) {
            // For each query base
            
            if(naturalMappings[i].isMapped()) {
                // If it actually mapped...
                
                // Work out where to
                TextPosition candidate = naturalMappings[i].getLocation();
                
                // Make a merge to that position.
                generateMerge(queryContig, i + 1, 
                    candidate.getContigNumber(), index.getOffset(candidate),
                    candidate.getStrand());
                    
                mappedBases++;
            }
        }
        
        if(credit) {
            // Each base zips matching bases inwards until it hits a mismatch,
            // or another mapped base. If zippings disagree, the base is not to
            // be mapped. Only bases between pairs of mapped bases can be
            // mapped.
            
            Log::info() << "Applying credit tolerating " << z_max <<
                " mismatches." << std::endl;
            
            // Make a set of all the zippings we will find for each position.
            std::vector<std::set<TextPosition>> zippings(
                naturalMappings.size());
        
            // Find the first mapped base (past the end if none)
            size_t firstMapped = 0;
            for(; firstMapped < naturalMappings.size() && 
                !naturalMappings[firstMapped].isMapped(); firstMapped++);
            
            // Find the last mapped base (size_t version of -1 if none)
            size_t lastMapped = naturalMappings.size() - 1;
            for(; lastMapped != (size_t) -1 && 
                !naturalMappings[lastMapped].isMapped(); lastMapped--);
            
            // This holds the index of the base currently providing credit.
            size_t provider;
            
            // How many mismatches have we seen since the credit provider?
            size_t mismatchesSeen;
            
            for(size_t i = firstMapped; i < naturalMappings.size() &&
                i <= lastMapped; i++) {
                
                // From the first to the last, do all the zipping off the right
                // sides of mapped bases.
                
                if(naturalMappings[i].isMapped()) {
                    // This base is mapped and is now the credit provider
                    provider = i;
                    mismatchesSeen = 0;
                } else {
                    // This base is not mapped. We know we saw a mapped base
                    // already and thus have a credit provider.
                    
                    // What position would we be implied to be at?
                    TextPosition implied = 
                        naturalMappings[provider].getLocation();
                    implied.addLocalOffset((int64_t) i - (int64_t) provider);
                    
                    if(implied.getOffset() >= index.getContigLength(
                        implied.getContigNumber())) {
                        
                        // This position is implied to zip off the end of the
                        // reference text. Skip to the next one.
                        continue;
                    }
                        
                    if(getCharacter(index, implied) != contig[i]) {
                        // This is a mismatch
                        mismatchesSeen++;
                        
                        Log::trace() << getCharacter(index, implied) << 
                            " at " << implied << " mismatches " << contig[i] <<
                            std::endl;
                    }
                    
                    if(mismatchesSeen <= z_max) {
                        // Not too many mismatches since the last credit
                        // provider. We can do credit.
                        
                        // TODO: not checking base identity here lets credit
                        // conflict even if it's wrong about base identity.
                        
                        // Zip this base to the position it is implied as being
                        // at by the credit provider.
                        zippings[i].insert(implied);
                        
                        Log::trace() << "Right credit zips " << i << " to " <<
                            implied << std::endl;
                        
                    }    
                }
                
            }
            
            for(size_t i = lastMapped; i != (size_t) -1 && i >= firstMapped;
                i--) {
                
                // From the last to the first, do all the zipping off the left
                // sides of mapped bases.
                
                // TODO: Unify this stateful logic with the above; it's
                // direction-independent.
                
                if(naturalMappings[i].isMapped()) {
                    // This base is mapped and is now the credit provider
                    provider = i;
                    mismatchesSeen = 0;
                } else {
                    // This base is not mapped. We know we saw a mapped base
                    // already and thus have a credit provider.
                    
                    // What position would we be implied to be at?
                    TextPosition implied = 
                        naturalMappings[provider].getLocation();
                    implied.addLocalOffset((int64_t) i - (int64_t) provider);
                    
                    if(implied.getOffset() >= index.getContigLength(
                        implied.getContigNumber())) {
                        
                        // This position is implied to zip off the end of the
                        // reference text. Skip to the next one.
                        continue;
                    }
                        
                    if(getCharacter(index, implied) != contig[i]) {
                        // This is a mismatch
                        mismatchesSeen++;
                        
                        Log::trace() << getCharacter(index, implied) <<
                            " at " << implied << " mismatches " << contig[i] <<
                            std::endl;
                    }
                    
                    if(mismatchesSeen <= z_max) {
                        // Not too many mismatches since the last credit
                        // provider. We can do credit.
                        
                        // TODO: not checking base identity here lets credit
                        // conflict even if it's wrong about base identity.
                        
                        // Zip this base to the position it is implied as being
                        // at by the credit provider.
                        zippings[i].insert(implied);
                        
                        Log::trace() << "Left credit zips " << i << " to " <<
                            implied << std::endl;
                    }    
                }
                
                // We don't need to do another pass, we can integrate here.
                
                if(zippings[i].size() == 1) {
                    // This base got zipped to exactly one place.
                    
                    // Work out where to (the only element in the set)
                    TextPosition candidate = *(zippings[i].begin());
                    
                    if(getCharacter(index, candidate) == contig[i]) {
                        // Make a merge to that position.
                        generateMerge(queryContig, i + 1, 
                            candidate.getContigNumber(),
                            index.getOffset(candidate), candidate.getStrand());
                            
                        Log::debug() << "Credit agrees on " << i << std::endl;
                        
                        mappedBases++;
                        creditBases++;     
                    } else {
                        // Merging to that position would merge mismatching
                        // bases.
                        Log::debug() << "Credit agrees on " << i <<
                            " but would merge a mismatch." << std::endl;
                    }
                        
                    
                        
                } else if(zippings[i].size() > 1) {
                    // This base has conflicted credit.
                    
                    Log::debug() << "Credit disagrees on " << i << std::endl;
                    
                    conflictedCredit++;
                }
                
            }
        }
        
        Log::output() << "Mapped " << mappedBases << " bases, " << 
            creditBases << " on credit, " << conflictedCredit << 
            " bases with conflicting credit." << std::endl;
    }
}







