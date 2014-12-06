#include "LRMappingScheme.hpp"
#include "Log.hpp"
#include "CreditFilter2.hpp"
#include "DisambiguateFilter.hpp"

#include <vector>

void LRMappingScheme::map(const std::string& query,
    std::function<void(size_t, TextPosition)> callback) const {
    
    // How many bases have we mapped or not mapped (credit or not).
    size_t mappedBases = 0;
    size_t unmappedBases = 0;
    
    // How many bases are mapped on credit?
    size_t creditBases = 0;
    // How many bases have conflicted credit?
    size_t conflictedCredit = 0;
    
    // How many are mapped on both the left and the right?
    size_t leftRightMappedBases = 0;
        
    // We will fill these in with mappings, which automatically come out
    // populated with TextPositions that are in the appropriate ranges.
    std::vector<Mapping> rightMappings;
    std::vector<Mapping> leftMappings;    
    
    // Map it on the right.
    rightMappings = index.misMatchMap(ranges, query, mask, minContext,
        addContext, multContext, z_max);
    
    if(addContext == 0 && multContext == 0 && z_max == 0) {
     
        // Take a moment to make sure we are working properly wrt the exact
        // match mappings.
        
        std::vector<Mapping> otherRightMappings = index.mapRight(query,
            mask, minContext);
            
        for(size_t i = 0; i < rightMappings.size(); i++) {
            if(rightMappings[i] != otherRightMappings[i]) {
                throw std::runtime_error("Mapping mismatch!");
            }
        }
        
    }
    
    // Map it on the left
    leftMappings = index.misMatchMap(ranges, reverseComplement(query), mask,
        minContext, addContext, multContext, z_max); 
    
    // Flip the left mappings back into the original order. They should stay
    // as other-side ranges.
    std::reverse(leftMappings.begin(), leftMappings.end());
    
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
        filteredMappings = CreditFilter2(index, ranges, z_max, 
            mask).apply(leftMappings, rightMappings,
            query);
    } else {
        // Apply a disambiguate filter to the mappings
        filteredMappings = DisambiguateFilter(index).apply(leftMappings,
            rightMappings);
    }
    
    Log::output() << "LR Mappings:" << std::endl;
    
    for(size_t i = 0; i < query.size(); i++ ) {
        Log::output() << query[i] << "\t" << leftMappings[i] << "\t" << 
            filteredMappings[i] << "\t" << rightMappings[i] << std::endl;
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
        char queryBase = query[i];
        
        // Make sure we have the contig we are mapping to.
        const std::string& referenceContig = index.displayContigCached(
            candidate.getContigNumber());
        
        // Grab the base from the contig.
        // TODO: factor this into FMDIndex as like displayCached or something.
        char referenceBase;
        if(candidate.getStrand()) {
            // Flip around a copy of the position and pull the complement of
            // what's on the forward strand, so we don't go and complement the
            // whole contig.
            TextPosition flipped = candidate;
            flipped.flip(index.getContigLength(candidate.getContigNumber()));
            
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
        
            // Dispatch the callback with this query index and this
            // TextPosition.
            callback(i, candidate);
            
            // Note we mapped a base
            mappedBases++;
            
            if(leftMappings[i].isMapped() && rightMappings[i].isMapped()) {
                // This base was left-mapped and right-mapped (and also not
                // mapped on credit).
                leftRightMappedBases++;
            }
            
            // TODO: provide a way to dump Mapping objects.
            
        } else {
            // This base didn't map
            unmappedBases++;
        }
        
    }
    
}
