#include "MarkovModel.hpp"

#include "Log.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include <boost/algorithm/string.hpp>

MarkovModel::MarkovModel(std::string filename): logProbabilities(), order() {

    // We need to load up the kmers file, store all the counts, and calculate
    // all the log probabilities for the last characters.
    
    // We'll fill in this map with the parsed counts. In general the counts
    // might not be ints. This maps from prefix to a map from next character to
    // count.
    std::map<std::string, std::map<char, double>> kmerCounts;
    
    // Open the file for reading.
    std::ifstream file(filename);
    
    // We're going to break the file into lines
    std::string line;
    while(std::getline(file, line)) {
        // For each line, split it into parts
        std::vector<std::string> parts;
        
        // Split on spaces of any kind
        boost::split(parts, line, boost::is_space());
        
        if(parts.size() == 1 && parts[0] == "") {
            // Skip blank lines
            continue;
        }
        
        if(parts.size() != 2) {
            // Complain we got a bad input file.
            throw std::runtime_error(std::string("Invalid number of parts: ") +
                line);
        }
        
        if(parts[0].size() < 1) {
            // Don't take empty kmers. TODO: make sure they all have constant
            // length.
            throw std::runtime_error("Got a too-short kmer!");
        }
        
        if(kmerCounts.size() == 0) {
            // This is our very first item. Autodetect the order.
            order = parts[0].size() - 1;
        } else {
            if(parts[0].size() - 1 != order) {
                // Complain that the model doesn't know what order it is.
                throw std::runtime_error("Model order is inconsistent");
            }
        }
        
        // Grab the prefix from the kmer (possibly "")
        std::string prefix = parts[0].substr(0, parts[0].size() - 1);
        
        // And the character that comes after it.
        char nextChar = parts[0][parts[0].size() - 1];
        
        if(!kmerCounts.count(prefix)) {
            // We need to make a map for this prefix
            kmerCounts[prefix] = std::map<char, double>();
        }
        
        // Parse out the count and save it.
        kmerCounts[prefix][nextChar] = strtod(parts[1].c_str(), NULL);
        
        
    }
    
    // OK now we loaded the kmer counts, do the by-prefix normalization.
    for(auto prefixPair : kmerCounts) {
        Log::debug() << "Normalizing prefix " << prefixPair.first << std::endl;
        
        // We'll sum up the counts.
        double total = 0;
        
        for(auto nextCharPair : prefixPair.second) {
            // Sum up the total probability
            total += nextCharPair.second;
        }
        
        
        for(auto nextCharPair : prefixPair.second) {
            // Fill in the log probabilities
            double logProbability = std::log2(nextCharPair.second / total);
            logProbabilities[prefixPair.first + nextCharPair.first] = 
                logProbability;
        }
        
    }
    
}

double MarkovModel::encodingCost(const std::string& prefix, char next) {
    if(prefix.size() < order) {
        // We would need to go into start cahracters here, and that doesn't
        // actually make sense for the way we want to use this model (i.e.
        // sequences/reads can start in the middle of what we trained on). TODO:
        // train on a bunch of reads and/or add pseudocounts so all starts are
        // possible.
        
        // For now just say we don't cost anything until we can really check.
        return 0;
    }
    
    // Pull off the right number of characters from the end of the string, tack
    // on the next character, and work out how probable this combination is.
    auto found = logProbabilities.find(prefix.substr(prefix.size() - order, 
        order) + next);
    
    if(found == logProbabilities.end()) {
        // We never saw this at all.
        // Log probability would be -inf for impossible, it costs inf to encode.
        return std::numeric_limits<double>::infinity();
    } else {
        // Return the negation of what we found (since the log probs are
        // negative bits, and the encoding costs are positive bits).
        return -(found->second);
    }
}

double MarkovModel::encodingCost(const std::string& subtext) {
    // We're going to total up the cost here.
    double total = 0;
    
    for(size_t i = order; i < subtext.size(); i++) {
        // Sum up the encoding cost for each character in turn. We don't look
        // before we have at least order characters to use.
        total += encodingCost(subtext.substr(i - order, order), subtext[i]);
    }
    
    return total;
}

// What start/stop character is used?
const std::string MarkovModel::START_STOP = "=";

















