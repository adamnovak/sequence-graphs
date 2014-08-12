#ifndef FASTA_HPP
#define FASTA_HPP

#include <string>
#include <iostream>
#include <istream>
#include <fstream>

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
     * Returns the next FASTAS record. May not be called if getNext is false.
     */
    std::string getNext();
    
protected:
    // Keeps track of the open file.
    std::ifstream stream;
    
};

#endif
