#ifndef FASTA_HPP
#define FASTA_HPP

#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <utility>

/**
 * Defines a FASTA file which can be read from disk and iterated over (in a
 * nonstandard way because files are hard).
 */
class Fasta {
public:
    
    /**
     * Make a new FASTA by opening the given filename.
     */
    Fasta(std::string filename);
    
    /**
     * Returns true if there is a next FASTA record, false otherwise.
     */
    bool hasNext();
    
    /**
     * Returns the next FASTA sequence. May not be called if hasNext is false.
     */
    std::string getNext();
    
    /**
     * Returns the next FASTA record (ID and sequence). May not be called if
     * hasNext is false. Discards the description, if one is present.
     */
    std::pair<std::string, std::string> getNextRecord();
    
protected:
    // Keeps track of the open file.
    std::ifstream stream;
    
};

#endif
