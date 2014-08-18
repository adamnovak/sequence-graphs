#include <stdexcept>

#include <Log.hpp>
#include <Util.h> // From libsuffixtools, for reverse_complement

#include "MappingMergeScheme.hpp"

const size_t MappingMergeScheme::MAX_THREADS = 32;

MappingMergeScheme::MappingMergeScheme(const FMDIndex& index, 
    const BitVector& rangeVector, 
    const std::vector<std::pair<std::pair<size_t, size_t>, bool> >& rangeBases, 
    const BitVector& includedPositions, size_t genome, size_t minContext,
    bool credit, bool mismatch, size_t z_max) : 
    MergeScheme(index), threads(), queue(NULL), contigsToMerge(NULL),
    rangeVector(rangeVector), rangeBases(rangeBases),
    includedPositions(includedPositions), genome(genome),
    minContext(minContext), credit(credit), mismatch(mismatch), z_max(z_max) {
    
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
        //threads.push_back(std::thread(&MappingMergeScheme::generateMerges,
        //        this, contigsToMerge));
        
        // TODO: for now, make all the merges at once.
        generateMerges(contigsToMerge);
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
	size_t creditBases = 0;
	size_t conflictedBases = 0;
        
        // Grab the contig as a string
        std::string contig = index.displayContig(queryContig);
        
        // How many positions are available to map to?
        Log::info() << taskName << " mapping " << contig.size() << 
            " bases via " << BitVectorIterator(includedPositions).rank(
            includedPositions.getSize()) << " bottom-level positions" <<
            std::endl;
	    
	std::vector<std::pair<int64_t,size_t>> rightMappings;
	std::vector<std::pair<int64_t,size_t>> leftMappings;    
	
	if (mismatch) {
      
	    Log::info() << "Using mismatch mapping with z_max " << z_max << std::endl;
	
	    // Map it on the right
	    rightMappings = index.misMatchMapRight(rangeVector, contig,
		&includedPositions, minContext, z_max);
	
	    // Map it on the left
	    leftMappings = index.misMatchMapRight(rangeVector, 
		reverseComplement(contig), &includedPositions, minContext, z_max); 
		
	} else {
		
	    // Map it on the right
	    rightMappings = index.mapRight(rangeVector, contig,
		    &includedPositions, minContext);
    
	    // Map it on the left
	    leftMappings = index.mapRight(rangeVector, reverseComplement(contig),
		    &includedPositions, minContext);
    
	}
        
        // Flip the left mappings back into the original order. They should stay
        // as other-side ranges.
        std::reverse(leftMappings.begin(), leftMappings.end());
          
	int64_t leftSentinel;
	int64_t rightSentinel;
	std::vector<size_t> creditCandidates;
	
	for(size_t i = 0; i < leftMappings.size(); i++) {
	    // Scan for the first left-mapped position in this contig
	    
	    if(leftMappings[i].first != -1) {
		if(rightMappings[i].first == -1) {
		    leftSentinel = i;
		    break;
		    
		} else if(rightMappings[i].first != -1 &&
		    rangeBases[leftMappings[i].first].first == rangeBases[rightMappings[i].first].first &&
		    rangeBases[leftMappings[i].first].second != rangeBases[rightMappings[i].first].second) {
		    
		    leftSentinel = i;
		    break;
		
		}
	    }
	}
	
	for(size_t i = rightMappings.size() - 1; i > 0; i--) {
	    // Scan for the last right-mapped position in this contig
	    
	    if(rightMappings[i].first != -1) {
		if(leftMappings[i].first == -1) {
		    rightSentinel = i;
		    break;
		    
		} else if(leftMappings[i].first != -1 &&
		    rangeBases[leftMappings[i].first].first == rangeBases[rightMappings[i].first].first &&
		    rangeBases[leftMappings[i].first].second != rangeBases[rightMappings[i].first].second) {
		    
		    rightSentinel = i;
		    break;
		
		}
	    }
	}
	
	Log::info() << "Left sentinel is " << leftSentinel << ", right sentinel is " << rightSentinel << std::endl;
	
	// Now identify and merged mapped bases from the individual left and
	// right mappings
	
	for(size_t i = 0; i < leftMappings.size(); i++) {
	    // For each position, look at the mappings.
	    
	    if(leftMappings[i].first != -1) {
		// We have a left mapping. Grab its base.
		auto leftBase = rangeBases[leftMappings[i].first];
		
		if(rightMappings[i].first != -1) {
		    // We have a right mapping. Grab its base too.
		    auto rightBase = rangeBases[rightMappings[i].first];
		    
		    // Compare the position (contig, base) pairs (first) and the
		    // orientation flags (second)
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
			Log::info() << "Anchor Merged pos " << i << ", a(n) " << contig[i] << " on contig " << queryContig << " to " << leftBase.first.second << " on contig " << leftBase.first.first << " with orientation " << !leftBase.second << std::endl;
			
			mappedBases++;                   
			} else {
			    // Didn't map this one
			    unmappedBases++;
			    conflictedBases++;
			    if(i > leftSentinel && rightSentinel > i) {
				creditCandidates.push_back(i);
			    }
			    Log::info() << "Conflicted " << i << " " << contig[i] << std::endl;
			}
			
		} else {
		    // Left mapped and right didn't.
		    
		    // Do the same thing, taking the left base and merging into it.
		    // Orientation is backwards to start with from our backwards
		    // right-semantics, so flip it.
		    generateMerge(queryContig, i + 1, leftBase.first.first, 
				  leftBase.first.second, !leftBase.second);
		    Log::info() << "Anchor Merged pos " << i << ", a(n) " << contig[i] << " on contig " << queryContig << " to " << leftBase.first.second << " on contig " << leftBase.first.first << " with orientation " << !leftBase.second << std::endl;
		    
		    
		    mappedBases++;
		}
		
	    } else if(rightMappings[i].first != -1) {
		// Right mapped and left didn't.
		
		// We have a right mapping. Grab its base too.
		auto rightBase = rangeBases[rightMappings[i].first];
		
		// Merge with the same contig and base. Leave the orientation alone
		// (since it's backwards to start with).
		generateMerge(queryContig, i + 1, rightBase.first.first, 
			      rightBase.first.second, rightBase.second);
		Log::info() << "Anchor Merged pos " << i << ", a(n) " << contig[i] << " on contig " << queryContig << " to " << rightBase.first.second << " on contig " << rightBase.first.first << " with orientation " << rightBase.second << std::endl;
		
		mappedBases++; 
	    } else {
		// Didn't map this one
		unmappedBases++;
		if(i > leftSentinel && rightSentinel > i) {
		    creditCandidates.push_back(i);
		    
		}
	    }   
	}
	
	if(credit) {
	    // Iterate across all possible positions for mapping on credit
	    
	    int64_t firstR;
	    bool contextMappedR;
	    int64_t firstL;
	    bool contextMappedL;
	    std::pair<std::pair<size_t,size_t>,bool> firstBaseR;
	    std::pair<std::pair<size_t,size_t>,bool> firstBaseL;    
	    int64_t LROffset;
	    
	    size_t maxRContext = 0;
	    size_t maxLContext = 0;
	    
	    Log::info() << "Checking " << creditCandidates.size() << " unmapped positions \"in the middle\"" << std::endl;
	    
	    // Scan for the maximum 
	    
	    for(size_t i; i < rightMappings.size(); i++) {
		if(rightMappings[i].second > maxRContext) {
		    maxRContext = rightMappings[i].second;
		}
		if(leftMappings[i].second > maxLContext) {
		    maxLContext = leftMappings[i].second;
		}
	    }
	    
	    for(size_t i = 0; i < creditCandidates.size(); i++) {
		firstR = -1;
		firstL = -1;
		contextMappedR = true;
		contextMappedL = true;
		firstBaseL = std::make_pair(std::make_pair(0,0),0);
		firstBaseR = std::make_pair(std::make_pair(0,0),0);
		
		for(size_t j = creditCandidates[i] - 1; j > 0 && creditCandidates[i] - j < maxRContext; j--) {
		    // Search leftward from each position until you find a position
		    // whose context includes creditCandidates[i]
		    
		    if(firstR == -1) {
			if(rightMappings[j].second > creditCandidates[i] - j) {
			    firstBaseR = rangeBases[rightMappings[j].first];
			    firstR = j;
			}
			
			// Then continue to search left. If a position has creditCandidates[i]
			// in its right context we want to see if subsequent positions map to
			// the same place
			
			// Terminate the search at the first base you come to which no longer
			// has our base in its right context. It is not possible that a position
			// farther to the left has our base in its right context and continues
			// to map to the same place as before
			
		    } else {
			if(rightMappings[j].second > creditCandidates[i] - j) {
			    if(rangeBases[rightMappings[j].first].first.second != firstBaseR.first.second  + firstR - j ||
				rangeBases[rightMappings[j].first].second != firstBaseR.second ||
				rangeBases[rightMappings[j].first].first.first != firstBaseR.first.first) {
				contextMappedR = false;
				break;
			    }
			}
		    }
		}
		
		for(size_t j = creditCandidates[i] + 1; j < leftMappings.size() && j - creditCandidates[i] < maxLContext; j++) {
		    if(firstL == -1) {
			if(leftMappings[j].second > j - creditCandidates[i]) {
			    firstBaseL = rangeBases[leftMappings[j].first];
			    firstL = j;
			}
		    } else {
			if(leftMappings[j].second > j - creditCandidates[i]) {
			    if(rangeBases[leftMappings[j].first].first.second != firstBaseL.first.second + j - firstL ||
				rangeBases[leftMappings[j].first].second != firstBaseL.second ||
				rangeBases[leftMappings[j].first].first.first != firstBaseL.first.first) {
				contextMappedL = false;
			    break;
				}
			}
		    }
		}
		
		if(firstR != -1 && contextMappedR) {
		    if(firstBaseR.second) {
			LROffset = firstR - (int64_t)creditCandidates[i];
		    } else {
			LROffset = (int64_t)creditCandidates[i] - firstR;
		    }
		    int64_t temp = (int64_t)firstBaseR.first.second + LROffset;
		    firstBaseR.first.second = (size_t)temp;
		    
		    if(firstL != -1 && contextMappedL) {
			if(!firstBaseL.second) {
			    LROffset = firstL - (int64_t)creditCandidates[i];
			} else {
			    LROffset = (int64_t)creditCandidates[i] - firstL;
			}
			int64_t temp = (int64_t)firstBaseL.first.second + LROffset;
			firstBaseL.first.second = (size_t)temp;
			
			if(firstBaseR.first == firstBaseL.first &&
			    firstBaseR.second != firstBaseL.second) {
			    generateMerge(queryContig, creditCandidates[i] + 1, firstBaseL.first.first, 
					  firstBaseL.first.second, !firstBaseL.second);
			    //Log::info() << "LR Credit Merged pos " << creditCandidates[i] << ", a(n) " << contig[creditCandidates[i]] << " on contig " << queryContig << " to " << firstBaseR.first.second << " on contig " << firstBaseR.first.first << " with orientation " << firstBaseR.second << std::endl;
			    mappedBases++;
			unmappedBases--;
			creditBases++;
			
			    }
		    } else {
			generateMerge(queryContig, creditCandidates[i] + 1, firstBaseR.first.first, 
				      firstBaseR.first.second, firstBaseR.second);
			//Log::info() << "R Credit Merged pos " << creditCandidates[i] << ", a(n) " << contig[creditCandidates[i]] << " on contig " << queryContig << " to " << firstBaseR.first.second << " on contig " << firstBaseR.first.first << " with orientation " << firstBaseR.second << std::endl;
			mappedBases++;
			unmappedBases--;
			creditBases++;
			
		    }
		} else if (firstL != -1 && contextMappedL) {
		    //Log::info() << "Before " << firstBaseL.first.second << std::endl;
		    if(!firstBaseL.second) {
			LROffset = firstL - (int64_t)creditCandidates[i];
		    } else {
			LROffset = (int64_t)creditCandidates[i] - firstL;
		    }
		    int64_t temp = (int64_t)firstBaseL.first.second + LROffset;
		    firstBaseL.first.second = (size_t)temp;
		    //Log::info() << "After " << firstBaseL.first.second << ", moved by " << LROffset << std::endl;
		    generateMerge(queryContig, creditCandidates[i] + 1, firstBaseL.first.first, 
				  firstBaseL.first.second, !firstBaseL.second);
		    //Log::info() << "L Credit Merged pos " << creditCandidates[i] << ", a(n) " << contig[creditCandidates[i]] << " on contig " << queryContig << " to " << firstBaseL.first.second << " on contig " << firstBaseL.first.first << " with orientation " << firstBaseL.second << std::endl;
		    mappedBases++;
		    unmappedBases--;
		    creditBases++;
		}
	    }
	}
	
	// Report that we're done.
	Log::info() << taskName << " finished (" <<
		mappedBases << "|" << unmappedBases << ")" << std::endl;
	Log::info() << conflictedBases <<
		" unmapped on account of conflicting left and right context matches." << std::endl;
	
	if(credit) {
	    Log::info() << mappedBases - creditBases << " anchor mapped; " <<
	    creditBases << " extension mapped" <<  std::endl;
	}
            
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










