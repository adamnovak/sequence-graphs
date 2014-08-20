#include "LCPArray.hpp"

// We need SAElems, so we might as well explicitly include their file.
#include <STCommon.h>

#include "Log.hpp"

#include <iterator>
#include <iostream>
#include <fstream>

LCPArray::LCPArray(const SuffixArray& suffixArray, const ReadTable& strings): values(), 
    psvs(), nsvs() {

    if(suffixArray.getSize() == 0) {
        // Just have a 0-length LCP if there are absolutely no suffixes.
        return;    
    }
    
    // We shall use the algorithm of Kasai et al. 2001: "Linear-Time Longest-
    // Common-Prefix Computation in Suffx Arrays and Its Applications"
    // See <http://www.cs.iastate.edu/~cs548/references/linear_lcp.pdf>
    
    // Make the LCP big enough to hold everything.
    values.resize(suffixArray.getSize());
    
    // We need to iterate through suffix array positions in order by text
    // position. We make a map from SAElem to rank, exploiting the fact that
    // std::map is sorted. But we need to tell it a better sort order.
    std::map<SAElem, size_t, std::greater<SAElem>> ranks;
    
    Log::info() << "Calculating ranks" << std::endl;
    
    for(size_t i = 0; i < suffixArray.getSize(); i++) {
        // Fill in all the ranks
        ranks[suffixArray.get(i)] = i;
    }
    
    Log::info() << "Applying Kasai's Algorithm" << std::endl;
    
    // The algorithm keeps this state as it goes through suffixes in rank order.
    size_t height = 0;
    
    for(auto entry : ranks) {
        // Go through all the element, rank pairs. Entry is a pair of a SAElem
        // and its rank.
        
        // Grab the suffix as a local.
        SAElem currentSuffix = entry.first;
        
        Log::trace() << currentSuffix << " at rank " << entry.second <<
            std::endl;
        
        if(entry.second == 0) {
            // Skip the very first because it has no predecessor
            values[entry.second] = 0;
            continue;
        }
        
        // Get the suffix before this one
        SAElem prevSuffix = suffixArray.get(entry.second - 1);
        
        Log::debug() << "LCP of " << currentSuffix << " @ " << entry.second << 
            " and " << prevSuffix <<  " @ " << entry.second - 1 << " is: " <<
            std::endl;
        
        Log::trace() << getFromSuffix(currentSuffix, 0, strings) << 
            " vs. " << getFromSuffix(prevSuffix, 0, strings) << std::endl;
        
        while(getFromSuffix(currentSuffix, height, strings) == 
            getFromSuffix(prevSuffix, height, strings) && 
            getFromSuffix(prevSuffix, height, strings) != '$') {
            
            // While this suffix and the previous one match, keep scanning. We
            // automatically define end-of-text characters to be distinct for
            // each text, and since you can't have the same suffix of *the same
            // text* twice, we can just say they never match anything.

            // Advance the height by 1 since we matched.
            height++;
            
            Log::trace() << getFromSuffix(currentSuffix, height, strings) << 
                " vs. " << getFromSuffix(prevSuffix, height, strings) << 
                std::endl;
            
        }
        
        // Store the LCP value
        values[entry.second] = height;
        
        Log::debug() << "LCP[" << entry.second << "]=" << height << std::endl;
        
        if(height > 0) {
            // If height isn't 0, dial it back. Not really sure how that
            // balances out, but they proved this would work with math (despite
            // mixing up their variables), so it ought to work in practice.
            height--;
        }
        
    }
    
    // Now we should have calculated the whole LCP. Reclaim some memory.
    ranks.clear(); 
    
    Log::info() << "Scanning for PSVs" << std::endl;
    
    // OK, now we need to construct the PSV/NSV indexes. The easiest way is with
    // scans.
    for(size_t i = 0; i < values.size(); i++) {
        // For each suffix. TODO: iterators?
        
        // What's the position of the most recent smaller value? If none is
        // found, we'll just use 0.
        size_t psv = 0;
        
        for(size_t j = i - 1; j != (size_t) -1; j--) {
            // For each position before i in descending order (stopping on
            // underflow)...
            
            if(values[j] < values[i]) {
                // We found a smaller value.
                psv = j;
                break;
            }
        }
        
        // Stick in the PSV pointer
        psvs.push_back(psv);
    }
    
    // TODO: Just 1 scan?
    
    Log::info() << "Scanning for NSVs" << std::endl;
    
    // Now the same for the NSV
    for(size_t i = 0; i < values.size(); i++) {
        // For each suffix
        
        // What's the position of the next smaller value? If none is found,
        // we'll just use 1 past the end.
        size_t nsv = values.size();
        
        for(size_t j = i + 1; j < values.size(); j++) {
            // For each position after i in ascending order.
            
            if(values[j] < values[i]) {
                // We found a smaller value
                nsv = j;
                break;
            }            
        }
        
        // Stick in the NSV pointer
        nsvs.push_back(nsv);
        
        
    }
    
}

LCPArray::LCPArray(const std::string& filename) : values(), psvs(), nsvs() {
    // Make a binary input stream.
    std::ifstream file(filename, std::ifstream::binary);
    
    // This is going to hold how many items should be in each vector.
    size_t arrayLength;
    
    // Read the number of items in platform-native byte order.
    file.read((char*) &arrayLength, sizeof(size_t));
    
    // Resize all the vectors
    values.resize(arrayLength);
    psvs.resize(arrayLength);
    nsvs.resize(arrayLength);
    
    // Read in that many elements for each vector
    file.read((char*) &values[0], arrayLength * sizeof(size_t));
    file.read((char*) &psvs[0], arrayLength * sizeof(size_t));
    file.read((char*) &nsvs[0], arrayLength * sizeof(size_t));
    
    // Close up the file.
    file.close();
}

void LCPArray::save(const std::string& filename) const {
    // Make a binary output stream.
    std::ofstream file(filename, std::ios::out | std::ofstream::binary);
    
    // Grab the array length as a local
    size_t arrayLength = values.size();
    
    // Save the array length in platform-native byte order.
    file.write((char*) &arrayLength, sizeof(size_t));
    
    // Read in that many elements for each vector
    file.write((char*) &values[0], arrayLength * sizeof(size_t));
    file.write((char*) &psvs[0], arrayLength * sizeof(size_t));
    file.write((char*) &nsvs[0], arrayLength * sizeof(size_t));
    
    // Close up the file
    file.close();
    
}
