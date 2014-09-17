#ifndef MARKOVMODEL_HPP
#define MARKOVMODEL_HPP
/**
 * Markov.hpp: implementation of Markov model loading and scoring.
 */


#include <string>
#include <map>

class MarkovModel {

public:

    /**
     * Load a Markov model from a file. The file format is two whitespace-
     * separated columns, one of which contains kmers (of length 1 greater than
     * the order of the model) and the other of which contains counts. The
     * character '=' is used to represent the start and stop characters.
     *
     * All possible kmers must be included, including any desired pseudocounts.
     */
    MarkovModel(std::string filename);
    
    /**
     * Get the encoding cost for encoding the specified character at the end of
     * the given string.
     */
    double encodingCost(const std::string& prefix, char next);
    
    // What's the start/stop character?
    static const std::string START_STOP;
    
protected:
    
    // Holds log probability (negative, base 2, so measured in bits) for each
    // kmer. The probability of a kmer is defined as its count over the total
    // count of all kmers that match its first k-1 characters.
    std::map<std::string, double> logProbabilities;

};

#endif
