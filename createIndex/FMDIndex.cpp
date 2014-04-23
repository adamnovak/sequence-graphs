#include <iostream>
#include <cstdlib>
#include <numeric>
#include <sstream>

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

bool FMDIndex::getStrand(CSA::pair_type base) const {
    // What strand corresponds to that text? Strands are forward, then reverse.
    return base.first % 2 == 1;
}

CSA::usint FMDIndex::getOffset(CSA::pair_type base) const {
    // What base offset, 1-based, from the left, corresponds to this pair_type,
    // which may be on either strand.
    if(getStrand(base) == 0) {
        // We're on the forward strand. Make offset 1-based.
        return base.second + 1;
    } else {
        // We're on the reverse strand, so we measured from the end. Make it
        // 1-based.
        return lengths[getContigNumber(base)] - base.second;
    }
}

std::string FMDIndex::getName(CSA::pair_type base) const {

    // Unpack the coordinate parts.
    CSA::usint contig = getContigNumber(base);
    CSA::usint offset = getOffset(base);
    
    // Work out what to name the position.
    std::stringstream nameStream;
    // Leave off the strand.
    nameStream << "N" << contig << "B" << offset;
    return nameStream.str(); 
    
}

CSA::usint FMDIndex::getTotalLength() const {
    // Sum all the contig lengths and double (to make it be for both strands).
    // See <http://stackoverflow.com/a/3221813/402891>
    return std::accumulate(lengths.begin(), lengths.end(), 0) * 2;
}

char FMDIndex::display(CSA::usint contig, CSA::usint offset,
    bool strand) const {

    std::cout << "Displaying " << contig << ":" << offset << "." << strand << std::endl;

    // What offset into the text do we use?
    CSA::usint textOffset;
    
    if(!strand) {
        // On the forward strand. Keep the original offset, but convert to
        // 0-based.
        textOffset = offset - 1;
    } else {
        // On the reverse strand. Use the offset from the other end. Convert to
        // 0-based.
        textOffset = lengths[contig] - offset;
    }
    
    // After converting to a base in (text, offset) form, call the other
    // display.
    return display(std::make_pair(contig * 2 + strand, textOffset));
    
    
}

char FMDIndex::display(CSA::pair_type base) const {
    std::cout << "Displaying " << base.first << "," << base.second << std::endl;
    
    // Display 1 character off the appropriate strand.
    CSA::uchar* displayData = fmd.display(base.first, 
        std::make_pair(base.second, base.second));
        
    // Get the character value
    char toReturn = displayData[0];
    
    // Free the buffer that got allocated.
    free(displayData);
    
    // Send back the character.
    return toReturn;
}

