#include <iostream>
#include <cstdlib>
#include <numeric>

#include "FMDIndex.hpp"

FMDIndex::FMDIndex(std::string basename): fmd(basename), names(), lengths() {
    // We already loaded the index itself in the initializer. Go load the
    // length/order metadata.
    
    // Open the contig name/length file for reading.
    std::ifstream contigFile((basename + ".chrom.sizes").c_str());
    
    // Have a string to hold each line in turn.
    std::string line;
    while(std::getline(contigFile, line)) {
        // For each <contig>\t<length> pair...
        
        // Find the \t
        size_t separator = line.find('\t');
        
        // Split into name and length
        std::string contigName = line.substr(0, separator - 1);
        std::string contigLength = line.substr(separator + 1, line.size() - 1);
        
        // Parse the length
        long long int lengthNumber = atoll(contigLength.c_str());
        
        // Add it to the vector of names in number order.
        names.push_back(contigName);
        
        // And the vector of sizes in number order
        lengths.push_back(lengthNumber);
    }
    // Close up the contig file. We read our contig metadata.
    contigFile.close();
    
}

CSA::usint FMDIndex::getContigNumber(CSA::pair_type base) const {
    // What contig corresponds to that text? Contigs all have both strands.
    return base.first / 2;
}

CSA::usint FMDIndex::getStrand(CSA::pair_type base) const {
    // What strand corresponds to that text? Strands are forward, then reverse.
    return base.first % 2;
}

CSA::usint FMDIndex::getOffset(CSA::pair_type base) const {
    // What base offset, 1-based, from the left, corresponds to this pair_type,
    // which may be on either strand.
    if(getStrand(base) == 0) {
        // We're on the forward strand; just correct for the 1-basedness offset.
        return base.second + 1;
    } else {
        // We're on the reverse strand, so we measured from the end.
        return lengths[getContigNumber(base)] - base.second;
    }
}

CSA::usint FMDIndex::getTotalLength() const {
    // Sum all the contig lengths and double (for both strands). See
    // <http://stackoverflow.com/a/3221813/402891>
    return std::accumulate(lengths.begin(), lengths.end(), 0) * 2;
}

