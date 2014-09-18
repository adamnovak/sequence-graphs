#ifndef MARKOVMODEL_HPP
#define MARKOVMODEL_HPP
/**
 * Markov.hpp: implementation of Markov model loading and scoring.
 */


#include <string>
#include <map>


/**
 * Implements a Markov model.
 */
class MarkovModel {

protected:

    /**
     * Represents a Markov model state, and includes all the places we can go
     * next and the probabilities for each.
     */
    struct MarkovNode {
        // Where should we go on each character.
        // Most of this table will be unused but that is OK.
        MarkovNode* nextState[256] = {0,};
        
        // What's the log probability of each transition? The probability of a
        // kmer is defined as the count of the k+1-mer it corresponds to over
        // the total count of all kmers that match its first k-1 characters.
        double logProbability[256] = {0,};
    };

public:

    /**
     * You iterate through the model (in a thread-safe way) by keeping state
     * with one of these.
     */
    typedef MarkovNode* iterator;

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
    
    /**
     * Get the encoding cost for a subtext (which may not nexessarily begin
     * or end at sequence beginnign or ending points).
     */
    double encodingCost(const std::string& subtext);
    
    /**
     * Start in the Markov model at the end of the given history. If the history
     * is not sufficiently long to put you into a state in the model, returns
     * NULL.
     */
    iterator start(const std::string& history);
    
    /**
     * Retrun the encoding cost in bits for encoding the given character when in
     * the given state.
     */
    double encodingCost(iterator& state, char next);
    
    // What's the start/stop character?
    static const std::string START_STOP;
    
protected:

    // We keep our interlinked data structure of MarkovNodes in a map so we
    // don't lose any. The table pointers just let us go real fast through them.
    std::map<std::string, MarkovNode> nodes;
    
    // What order is the model (1 less than kmer length)
    size_t order;
    
private:
    // No copying or moving due to our using pointers we would have to do work
    // to update. TODO: implement move and move our node map without *actually*
    // changing its addresses internally. I am pretty sure that that is how move
    // works, but I haven't seen anything guaranteeing it.
    MarkovModel(const MarkovModel& other) = delete;
    MarkovModel(MarkovModel&& other) = delete;
    MarkovModel& operator=(const MarkovModel& other) = delete;
    MarkovModel& operator=(MarkovModel&& other) = delete;

};

#endif
