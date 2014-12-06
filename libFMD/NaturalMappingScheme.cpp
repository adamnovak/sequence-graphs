#include "NaturalMappingScheme.hpp"
#include "Log.hpp"

#include <vector>

void NaturalMappingScheme::map(const std::string& query,
    std::function<void(size_t, TextPosition)> callback) const {
    
    // Map using the natural context scheme: get matchings from all the
    // unique-in-the-reference strings that overlap you.
    
    // Call into the index. TODO: pass a parameters struct of some type.
    // TODO: It would make sense to use a DisambiguateFilter or some such,
    // but with this algorithm that's sort of integrated into the mapping.
    // TODO: This doesn't even use ranges.
    std::vector<Mapping> naturalMappings = index.naturalMap(query, mask,
        minContext);
        
    // How many bases have we mapped or not mapped (credit or not).
    size_t mappedBases = 0;
    size_t unmappedBases = 0;
    
    // How many bases are mapped on credit?
    size_t creditBases = 0;
    // How many bases have conflicted credit?
    size_t conflictedCredit = 0;
    
    for(size_t i = 0; i < naturalMappings.size(); i++) {
        // For each query base
        
        if(naturalMappings[i].isMapped()) {
            // If it actually mapped...
            
            // Work out where to
            TextPosition candidate = naturalMappings[i].getLocation();
            
            // Dispatch the callback with this query index and this
            // TextPosition.
            callback(i, candidate);
                
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
                    
                if(index.displayCached(implied) != query[i]) {
                    // This is a mismatch
                    mismatchesSeen++;
                    
                    Log::trace() << index.displayCached(implied) << 
                        " at " << implied << " mismatches " << query[i] <<
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
                    
                if(index.displayCached(implied) != query[i]) {
                    // This is a mismatch
                    mismatchesSeen++;
                    
                    Log::trace() << index.displayCached(implied) <<
                        " at " << implied << " mismatches " << query[i] <<
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
                
                if(index.displayCached(candidate) == query[i]) {
                    // Dispatch the callback with this query index and this
                    // TextPosition.
                    callback(i, candidate);
                        
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
