#include <stdexcept>

#include <Log.hpp>
#include <Util.h> // From libsuffixtools, for reverse_complement

#include "MappingMergeScheme.hpp"



MappingMergeScheme::MappingMergeScheme(const FMDIndex& index, 
    const BitVector& rangeVector, 
    const std::vector<std::pair<std::pair<size_t, size_t>, bool> >& rangeBases, 
    const BitVector& includedPositions, size_t genome, size_t minContext, bool credit, std::string mapType) : 
    MergeScheme(index), threads(), queue(NULL), rangeVector(rangeVector), 
    rangeBases(rangeBases), includedPositions(includedPositions), 
    genome(genome), minContext(minContext), credit(credit), mapType(mapType) {
    
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
	if(mapType == "LRexact") {
	    threads.push_back(std::thread(&MappingMergeScheme::generateMerges,
		    this, contig));
	    
	} else if(mapType == "centered") {
	    threads.push_back(std::thread(&MappingMergeScheme::CgenerateMerges,
		    this, contig));
	    
	} else {
	    throw std::runtime_error(mapType + " is not an implemented context type.");
	  
	}
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

void MappingMergeScheme::CgenerateMerges(size_t queryContig) const {
    
    // What's our thread name?
    std::string threadName = "T" + std::to_string(genome) + "." + 
        std::to_string(queryContig);
        
    // How many bases have we mapped or not mapped
    size_t mappedBases = 0;
    size_t unmappedBases = 0;
    size_t creditBases = 0;
    
    // Grab the contig as a string
    std::string contig = index.displayContig(queryContig);
    
    // How many positions are available to map to?
    Log::info() << threadName << " mapping " << contig.size() << 
        " bases via " << BitVectorIterator(includedPositions).rank(
        includedPositions.getSize()) << " bottom-level positions" << std::endl;
    
    // Map it
    std::vector<std::pair<int64_t,size_t>> Mappings = index.Cmap(rangeVector, contig, 
        &includedPositions, minContext);
    
    size_t leftSentinel;
    size_t rightSentinel;
    std::vector<size_t> creditCandidates;
    
    for(size_t i = 0; i < Mappings.size(); i++) {
	// Scan for the first mapped position in this contig
      
	if(Mappings[i].first != -1) {
	    leftSentinel = i;
	    break;

	}
    }
            
    for(size_t i = Mappings.size() - 1; i > 0; i--) {
	// Scan for the last mapped position in this contig

	if(Mappings[i].first != -1) {
	    leftSentinel = i;
	    break;

	}
    }      
       
      
    for(size_t i = 0; i < Mappings.size(); i++) {
        // For each position, look at the mappings.
        
        if(Mappings[i].first != -1) {
            // We have a mapping. Grab its base.
            auto Base = rangeBases[Mappings[i].first];
            generateMerge(queryContig, i + 1, Base.first.first, 
                        Base.first.second, !Base.second); 
                        
            mappedBases++;
	    Log::info() << "Mapped " << i << ", base " << contig[i] << std::endl;
	    
	} else {
            // Didn't map this one
            unmappedBases++;
	    creditCandidates.push_back(i);
	  
	}
    }
    
    bool firstR;
    bool contextMappedR;
    bool firstL;
    bool contextMappedL;
    std::pair<std::pair<size_t,size_t>,bool> firstBaseR;
    std::pair<std::pair<size_t,size_t>,bool> firstBaseL;    
    
    size_t maxContext = 0;
    
    Log::info() << "Checking " << creditCandidates.size() << " unmapped positions \"in the middle\"" << std::endl;
    
    for(size_t i; i < Mappings.size(); i++) {
	if(Mappings[i].second > maxContext) {
	    maxContext = Mappings[i].second;
	}
    }
        
    for(size_t i = 0; i < creditCandidates.size(); i++) {
	firstR = false;
	firstL = false;
	contextMappedR = true;
	contextMappedL = true;	
	for(size_t j = creditCandidates[i] - 1; j + 1 > 0 && creditCandidates[i] - j < maxContext; j--) {
	    // Search leftward from each position until you find a position
	    // whose context includes creditCandidates[i]
	  
	    if(!firstR) {
		if(Mappings[j].second > creditCandidates[i] - j) {
		    firstBaseR = rangeBases[Mappings[j].first];
		    firstR = true;
		}
		
	    // Then continue to search left. If a position has creditCandidates[i]
	    // in its right context we want to see if subsequent positions map to
	    // the same place
	    
	    // Terminate the search at the first base you come to which no longer
	    // has our base in its right context. It is not possible that a position
	    // farther to the left has our base in its right context and continues
	    // to map to the same place as before
		
	    } else {
		if(Mappings[j].second > creditCandidates[i] - j) {
		    if(rangeBases[Mappings[j].first].first.second != firstBaseR.first.second ||
			    rangeBases[Mappings[j].first].second != firstBaseR.second ||
			    rangeBases[Mappings[j].first].first.first != firstBaseR.first.first + creditCandidates[i] - j) {
			contextMappedR = false;
			break;
		    }
		}
	    }
	}
	
	for(size_t j = creditCandidates[i] + 1; j < creditCandidates.size() && j - creditCandidates[i] < maxContext; j++) {
	    if(!firstL) {
		if(Mappings[j].second > j - creditCandidates[i]) {
		    firstBaseL = rangeBases[Mappings[j].first];
		    firstL = true;
		}
	    } else {
		if(Mappings[j].second > j - creditCandidates[i]) {
		    if(rangeBases[Mappings[j].first].first.second != firstBaseL.first.second ||
			    rangeBases[Mappings[j].first].second != firstBaseL.second ||
			    rangeBases[Mappings[j].first].first.first != firstBaseL.first.first + j - creditCandidates[i]) {
			contextMappedL = false;
			break;
		    }
		}
	    }
	}
		    
		    
	if(firstR && contextMappedR) {
	    firstBaseR.first.first++;

	    if(firstL && contextMappedL) {
		firstBaseL.first.first--;
		
		if(firstBaseR.first == firstBaseL.first &&
			  firstBaseR.second != firstBaseL.second) {
		  generateMerge(queryContig, creditCandidates[i], firstBaseR.first.first, 
                        firstBaseR.first.second, firstBaseR.second);
		  Log::info() << "Credit Merged " << creditCandidates[i] << std::endl;
		  mappedBases++;
		  unmappedBases--;
		  creditBases++;
		  
		}
	    } else {
		  generateMerge(queryContig, creditCandidates[i], firstBaseR.first.first, 
                        firstBaseR.first.second, firstBaseR.second);
		  Log::info() << "Credit Merged " << creditCandidates[i] << std::endl;
		  mappedBases++;
		  unmappedBases--;
		  creditBases++;
		  
	    }
	} else if (firstL && contextMappedL) {
	    firstBaseL.first.first--;
	    generateMerge(queryContig, creditCandidates[i], firstBaseL.first.first, 
                        firstBaseL.first.second, !firstBaseL.second);
	    mappedBases++;
	    unmappedBases--;
	    creditBases++;
	    Log::info() << "Credit Merged " << creditCandidates[i] << std::endl;
	  
	}
		
    }
    
    // Close the queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
    // Report that we're done.
    Log::info() << threadName << " Anchor Merged (" << mappedBases << "|" << 
        unmappedBases << ")." << " Credit Merged " << creditBases << std::endl;
	
}

void MappingMergeScheme::generateMerges(size_t queryContig) const {
    
    // What's our thread name?
    std::string threadName = "T" + std::to_string(genome) + "." + 
        std::to_string(queryContig);
        
    // How many bases have we mapped or not mapped
    size_t mappedBases = 0;
    size_t unmappedBases = 0;
    size_t creditBases = 0;
    
    // Grab the contig as a string
    std::string contig = index.displayContig(queryContig);
    
    // How many positions are available to map to?
    Log::info() << threadName << " mapping " << contig.size() << 
        " bases via " << BitVectorIterator(includedPositions).rank(
        includedPositions.getSize()) << " bottom-level positions" << std::endl;
    
    // Map it on the right
    std::vector<std::pair<int64_t,size_t>> rightMappings = index.map(rangeVector, contig, 
        &includedPositions, minContext);
    
    // Map it on the left
    std::vector<std::pair<int64_t,size_t>> leftMappings = index.map(rangeVector, 
        reverseComplement(contig), &includedPositions, minContext);
    
    // Flip the left mappings back into the original order. They should stay as
    // other-side ranges.
    std::reverse(leftMappings.begin(), leftMappings.end());
    
     size_t leftSentinel;
    size_t rightSentinel;
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
                        
                    mappedBases++;                   
                } else {
                    // Didn't map this one
                    unmappedBases++;
		    if(i > leftSentinel && rightSentinel > i) {
			creditCandidates.push_back(i);
		    }
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
            
        } else if(rightMappings[i].first != -1) {
            // Right mapped and left didn't.
            
            // We have a right mapping. Grab its base too.
            auto rightBase = rangeBases[rightMappings[i].first];
            
            // Merge with the same contig and base. Leave the orientation alone
            // (since it's backwards to start with).
            generateMerge(queryContig, i + 1, rightBase.first.first, 
                        rightBase.first.second, rightBase.second); 
            
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
    
    bool firstR;
    bool contextMappedR;
    bool firstL;
    bool contextMappedL;
    std::pair<std::pair<size_t,size_t>,bool> firstBaseR;
    std::pair<std::pair<size_t,size_t>,bool> firstBaseL;    
    
    size_t maxRContext = 0;
    size_t maxLContext = 0;
    
    Log::info() << "Checking " << creditCandidates.size() << " unmapped positions \"in the middle\"" << std::endl;
    
    for(size_t i; i < rightMappings.size(); i++) {
	if(rightMappings[i].second > maxRContext) {
	    maxRContext = rightMappings[i].second;
	}
	if(leftMappings[i].second > maxLContext) {
	    maxLContext = leftMappings[i].second;
	}
    }
        
    for(size_t i = 0; i < creditCandidates.size(); i++) {
	firstR = false;
	firstL = false;
	contextMappedR = true;
	contextMappedL = true;	
	for(size_t j = creditCandidates[i] - 1; j + 1 > 0 && creditCandidates[i] - j < maxRContext; j--) {
	    // Search leftward from each position until you find a position
	    // whose context includes creditCandidates[i]
	  
	    if(!firstR) {
		if(rightMappings[j].second > creditCandidates[i] - j) {
		    firstBaseR = rangeBases[rightMappings[j].first];
		    firstR = true;
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
		    if(rangeBases[rightMappings[j].first].first.second != firstBaseR.first.second ||
			    rangeBases[rightMappings[j].first].second != firstBaseR.second ||
			    rangeBases[rightMappings[j].first].first.first != firstBaseR.first.first + creditCandidates[i] - j) {
			contextMappedR = false;
			break;
		    }
		}
	    }
	}
	
	for(size_t j = creditCandidates[i] + 1; j < creditCandidates.size() && j - creditCandidates[i] < maxLContext; j++) {
	    if(!firstL) {
		if(leftMappings[j].second > j - creditCandidates[i]) {
		    firstBaseL = rangeBases[rightMappings[j].first];
		    firstL = true;
		}
	    } else {
		if(leftMappings[j].second > j - creditCandidates[i]) {
		    if(rangeBases[leftMappings[j].first].first.second != firstBaseL.first.second ||
			    rangeBases[leftMappings[j].first].second != firstBaseL.second ||
			    rangeBases[leftMappings[j].first].first.first != firstBaseL.first.first + j - creditCandidates[i]) {
			contextMappedL = false;
			break;
		    }
		}
	    }
	}
		    
		    
	if(firstR && contextMappedR) {
	    firstBaseR.first.first++;

	    if(firstL && contextMappedL) {
		firstBaseL.first.first--;
		
		if(firstBaseR.first == firstBaseL.first &&
			  firstBaseR.second != firstBaseL.second) {
		  generateMerge(queryContig, creditCandidates[i], firstBaseR.first.first, 
                        firstBaseR.first.second, firstBaseR.second);
		  Log::info() << "Credit Merged " << creditCandidates[i] << std::endl;
		  mappedBases++;
		  unmappedBases--;
		  creditBases++;
		  
		}
	    } else {
		  generateMerge(queryContig, creditCandidates[i], firstBaseR.first.first, 
                        firstBaseR.first.second, firstBaseR.second);
		  Log::info() << "Credit Merged " << creditCandidates[i] << std::endl;
		  mappedBases++;
		  unmappedBases--;
		  creditBases++;
		  
	    }
	} else if (firstL && contextMappedL) {
	    firstBaseL.first.first--;
	    generateMerge(queryContig, creditCandidates[i], firstBaseL.first.first, 
                        firstBaseL.first.second, !firstBaseL.second);
	    mappedBases++;
	    unmappedBases--;
	    creditBases++;
	    Log::info() << "Credit Merged " << creditCandidates[i] << std::endl;
	  
	}
		
    }
    }
    
    // Close the queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
    // Report that we're done.
    Log::info() << threadName << " finished (" << mappedBases << "|" << 
        unmappedBases << ")" << std::endl;
	
    if(credit) {
	Log::info() << mappedBases - creditBases << " anchor mapped; " << creditBases << " extension mapped" <<  std::endl;
	}
}








