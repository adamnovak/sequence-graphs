#ifndef GENERICBITVECTOR_HPP
#define GENERICBITVECTOR_HPP

#include <string>
#include <fstream>
#include <istream>
#include <ostream>
#include <utility>
#include <mutex>

#define BITVECTOR_SDSL
#ifdef BITVECTOR_CSA
    // Using CSA bitvectors
    #include "BitVector.hpp"
#endif
#ifdef BITVECTOR_SDSL
    // Using SDSL bitvectors
    #include <sdsl/bit_vectors.hpp>
#endif

/**
 * Represents a bit vector that bastracts out its underlying implementation in
 * favor of a simple interface. Right now we're going to define it over RLCSA.
 *
 * Needs to support O(1) rank and select with a small constant factor.
 *
 * Needs to support multi-threaded rank/select access.
 */
class GenericBitVector {
public:
    /**
     * Load a bitvector from a stream.
     */
    GenericBitVector(std::ifstream& stream);

    /**
     * Load a bitvector from a stream we take.
     */
    GenericBitVector(std::ifstream&& stream);

    /**
     * Load a bitvector from a file.
     */
    GenericBitVector(const std::string& filename);
    
    /**
     * Create a bitvector from scratch. Must call finish before using rank and
     * select.
     */
    GenericBitVector();
    
    /**
     * Destroy a generic bitvector.
     */
    ~GenericBitVector();
    
    /**
     * Get the size of the bitvector in bits. Must have been finished first.
     */
    #ifdef BITVECTOR_CSA
    inline size_t getSize() const {
        return size;
    }
    #endif

    #ifdef BITVECTOR_SDSL
    inline size_t getSize() const {
        return bitvector.size();
    }
    #endif
    
    /**
     * Add a 1 bit at the given index. Must be called in increasing order.
     * Cannot be run on a bitvector that has been loaded or finished.
     */
    #ifdef BITVECTOR_CSA
    inline void addBit(size_t index) {
        if(encoder == NULL) {
            throw std::runtime_error("Can't add to a vector we didn't create!");
        }
        
        // Pass on valid requests.
        encoder->addBit(index);
    }
    #endif
    #ifdef BITVECTOR_SDSL
    inline void addBit(size_t index) {
        if(rankSupport != NULL || selectSupport != NULL) {
            throw std::runtime_error("Can't add to a finished/loaded vector!");
        }

        // What's the size of the bitvector before we make it bigger?
        size_t oldSize = getSize();

        if(index >= oldSize) {
            // Make sure we have room
            bitvector.resize(index + 1);
            
            for(size_t i = oldSize; i < index; i++) {
                // Make sure all the bits we just added in are 0s.
                bitvector[i] = 0;
            }
            
        }
        
        // Set the bit
        bitvector[index] = 1;
    }
    #endif
    
    /**
     * Returns true if the given index is set, or false otherwise. Must be
     * thread-safe.
     */
    #ifdef BITVECTOR_CSA
    inline bool isSet(size_t index) const {
        // Grab the iterator
        std::lock_guard<std::mutex> lock(iteratorMutex);
        
        return iterator->isSet(index);
    }
    #endif
    #ifdef BITVECTOR_SDSL
    inline bool isSet(size_t index) const {
        return bitvector[index];
    }
    #endif
    
    /**
     * Must be called on a bitvector not loaded from a file before using rank
     * and select, or saving to a file. Takes the total length of the bitvector
     * inclusing all trailing zeros.
     */
    void finish(size_t length);
    
    /**
     * Save the BitVector to the given stream. It hust have already been
     * finished.
     */
    void writeTo(std::ofstream& stream) const;
    
    /**
     * Get the number of 1s occurring before the given index. Must be thread-
     * safe.
     */
    #ifdef BITVECTOR_CSA
    inline size_t rank(size_t index) const {

        // Grab the iterator
        std::lock_guard<std::mutex> lock(iteratorMutex);
        
        // Take the rank of a position (not in at_least mode)
        return iterator->rank(index, false);
    }
    #endif
    #ifdef BITVECTOR_SDSL
    inline size_t rank(size_t index) const {
        if(index >= getSize()) {
            // Get the total number of 1s in the bitvector
            return (*rankSupport)(getSize()) + isSet(getSize() - 1);
        }
        
        // Take the rank of a position (not in at_least mode). Correct for
        // disagreement over what rank is.
        return (*rankSupport)(index + 1);
    }
    #endif

    
    /**
     * Get the number of 1s occurring before the given index. Must be thread-
     * safe. If at-least is set, includes a 1 at the given index.
     */
    inline size_t rank(size_t index, bool atLeast) const {
        if(atLeast) {
            // Implement this mode ourselves so we can be certain of exactly what it
            // does.
            if(index == 0) {
                // Do the calculation below, simplifying out the rank(-1) call which
                // would alwauys be 0 if we could actually run it.
                return (index < getSize());
            }
            // Get the rank of the previous position, then add 1 if this position
            // could contain a 1.
            return rank(index-1) + (index < getSize());
        } else {
            return rank(index);
        }
    }
    
    /**
     * Select the index at which the one with the given rank appears (i.e. the
     * first one would be 0). Must be thread-safe.
     */
    #ifdef BITVECTOR_CSA
    inline size_t select(size_t one) const {
        // Grab the iterator
        std::lock_guard<std::mutex> lock(iteratorMutex);
        
        // Go select the right position.
        return iterator->select(one);
    }
    #endif

    #ifdef BITVECTOR_SDSL
    inline size_t select(size_t one) const {
        // Go select the right position. We need to compensate since 0 means the
        // first one to us, but 1 means the first one to SDSL.
        return (*selectSupport)(one + 1);
    }
    #endif
    
    /**
     * OR two bitvectors together.
     */
    GenericBitVector* createUnion(const GenericBitVector& other) const;
    
    /**
     * Given an indes, return the index of the next 1 at or after that position,
     * paired with its rank.
     */
    #ifdef BITVECTOR_CSA
    inline std::pair<size_t, size_t> valueBefore(size_t index) const {
        // Get the rank of the 1 before this position, 1-based apparently.
        size_t lastRank = rank(index);
        if(lastRank == 0 || index >= getSize()) {
            // We ought to be giving the past-the-end value.
            lastRank = rank(getSize() + 1);
        } else {
            // Budge left to get the actual 0-based rank.
            lastRank--;
        }
        
        // Look up where it is and return that and the rank.
        return std::make_pair(select(lastRank), lastRank);
    }
    #endif
    #ifdef BITVECTOR_SDSL
    inline std::pair<size_t, size_t> valueBefore(size_t index) const {
        // Get the rank of the 1 before this position, 1-based apparently.
        size_t lastRank = rank(index);
        
        if(lastRank == 0) {
            // There is no value before, do the past-the-end result.
            return std::make_pair(getSize(), rank(getSize() - 1));
        }
        
        // Assume we are in bounds.
        
        // Budge left to get the actual 0-based rank.
        lastRank--;
        
        // Look up where it is and return that and the rank.
        return std::make_pair(select(lastRank), lastRank);
    }
    #endif
    
    /**
     * Given an indes, return the index of the last 1 at or before that
     * position, paired with its rank. Wraps around if no such value is found.
     */
    #ifdef BITVECTOR_CSA
    inline std::pair<size_t, size_t> valueAfter(size_t index) const {
        // Get the rank of the next 1 at or after the given position.
        size_t nextRank = rank(index, true) - 1;
        if(index >= getSize()) {
            // We ought to be giving the past-the-end value even though rank
            // dips by 1 at the end of the bitvector.
            nextRank += 1;
        }
        
        // Look up where it is and return that and the rank.
        return std::make_pair(select(nextRank), nextRank);
    }
    #endif
    #ifdef BITVECTOR_SDSL
    inline std::pair<size_t, size_t> valueAfter(size_t index) const {
        if(index >= getSize()) {
            // We ought to be giving the past-the-end value
            return std::make_pair(getSize(), rank(getSize() - 1));
        }
        // Assume we're in bounds
        
        // Get the rank of the next 1 at or after the given position.
        size_t nextRank = rank(index, true) - 1;
        
        if(nextRank >= rank(getSize() - 1)) {
            // Too close to the end to do a select. Return past the end value
            // again.
            return std::make_pair(getSize(), rank(getSize() - 1));
        }
        
        // Look up where it is and return that and the rank.
        return std::make_pair(select(nextRank), nextRank);
    }
    #endif
    
    // Move construction OK
    GenericBitVector(GenericBitVector&& other) = default;
    
    // Move assignment OK
    GenericBitVector& operator=(GenericBitVector&& other) = default; 
    
private:
    
    // No copy
    GenericBitVector(const GenericBitVector& other) = delete;
    
    // No assignment
    GenericBitVector& operator=(const GenericBitVector& other) = delete;
    
    #ifdef BITVECTOR_CSA
        // Actual implementation on CSA
        // We need an encoder to actually make the bitvector. We own this.
        BitVectorEncoder* encoder;
        // And we need a place to put the vector when done.
        BitVector* bitvector;
        // And we need to remember how far along we are in our encoding.
        size_t size;
        // And we need an iterator to query the BitVector when done.
        BitVectorIterator* iterator;
        // And a lock to protect that
        mutable std::mutex iteratorMutex;
    #endif
    #ifdef BITVECTOR_SDSL
        // Actual implementation on SDSL.
        
        // We seem to be only able to use full uncompressed bitvectors when 
        // we're actually building them...
        sdsl::bit_vector bitvector;
        
        // And we need rank on that
        sdsl::rank_support_v5<>* rankSupport;
        // And select
        sdsl::select_support_mcl<>* selectSupport;        
    #endif
    
};



#endif
