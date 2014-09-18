#include "MarkovModel.hpp"

#include "Log.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include <boost/algorithm/string.hpp>

MarkovModel::MarkovModel(std::string filename): nodes(), order() {

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
        
        // Make sure there is a node for each state that actually happens
        nodes[prefixPair.first] = MarkovNode();
        
        for(auto nextCharPair : prefixPair.second) {
            // Sum up the total probability
            total += nextCharPair.second;
            
        }
        
        
        for(auto nextCharPair : prefixPair.second) {
            // Work out the log probability for going to this node
            double logProbability = std::log2(nextCharPair.second / total);
            
            // And fill it in
            nodes[prefixPair.first].logProbability[nextCharPair.first] = 
                logProbability;
        }
        
    }
    
    for(auto prefixPair : kmerCounts) {
        // Now we need to fill in the pointers to the next states for each node.
        MarkovNode& node = nodes[prefixPair.first];
        
        for(auto nextCharPair : prefixPair.second) {
            // What does our memory look like when we add on this next
            // character?
            std::string nextStateName = prefixPair.first.substr(1, order) + 
                nextCharPair.first;
                
            // Make a pointer right there so we don't have to bother building
            // that string.
            node.nextState[nextCharPair.first] = &nodes[nextStateName];
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
    
    // What does our memory look like at the end of the string?
    std::string memory = prefix.substr(prefix.size() - order, order);
    
    // Jump to that node and get the log probability for this next character.
    // Negate it before returning.
    return -nodes[memory].logProbability[next];
}

double MarkovModel::encodingCost(const std::string& subtext) {
    if(order >= subtext.size()) {
        // No transitions observed. Give up now.
        return 0;
    }
    
    // We're going to total up the cost here.
    double total = 0;
    
    // Start now that we know we have enough state
    iterator state = start(subtext.substr(0, order));
    
    for(size_t i = order; i < subtext.size(); i++) {
        // Encode all the characters
        total += encodingCost(state, subtext[i]);
    }

    return total;
}

MarkovModel::iterator MarkovModel::start(const std::string& history) {
    if(history.size() < order) {
        // Can't get a state
        return NULL;
    } else if(history.size() == order) {
        // Just use this string for the lookup
        return &nodes[history];
    } else {
        // We need to pull off the end piece and use that
        std::string memory = history.substr(history.size() - order, order);
        return &nodes[memory];
    }
}

double MarkovModel::backfill(MarkovModel::iterator& state,
    const std::string& history) {
    
    if(order > history.size()) {
        // No transitions observed. Give up now.
        Log::debug() << "History " << history << " too short." << std::endl;
        state = NULL;
        return 0;
    }
    
    // We're going to total up the cost here.
    double total = 0;
    
    Log::debug() << "Backfilling from " << history.substr(0, order) <<
        std::endl;
    
    // Start now that we know we have enough state
    state = start(history.substr(0, order));
    
    for(size_t i = order; i < history.size(); i++) {
        // Encode all the characters, advancing the state
        double cost = encodingCost(state, history[i]);
        total += cost;
    }

    Log::debug() << "Backfilled " << total << " bits" << std::endl;

    // Give back the total, leaving the state where it should be.
    return total;
}

double MarkovModel::encodingCost(MarkovModel::iterator& state, char next) {
    // Remember where we are
    iterator oldState = state;
    
    // Go to the next state
    state = state->nextState[next];
    
    // Return the encoding cost to do so (negative log probability)
    double cost =-oldState->logProbability[next];
    Log::debug() << "Encoding " << next << " costs " << cost  << " bits" <<
        std::endl;
    return cost;
}

// What start/stop character is used?
const std::string MarkovModel::START_STOP = "=";

















