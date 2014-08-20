#ifndef GENERICBITVECTOR_HPP
#define GENERICBITVECTOR_HPP

#include <string>
#include <fstream>
#include <istream>
#include "BitVector.hpp"

/**
 * Represents a bit vector that bastracts out its underlying implementation in
 * favor of a simple interface. Right now we're going to define it over RLCSA.
 *
 * Needs to support O(1) rank and select with a small constant factor.
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
     * Must be called on a bitvector not loaded from a file before using rank
     * and select.
     */
    void finish();
    
    /**
     * Get the number of 1s occurring before the given index.
     */
    size_t rank(size_t index);
    
    /**
     * Select the index at which the one with the given rank appears (i.e. the
     * first one would be 0).
     */
    size_t select(size_t one);
    
    // Move construction OK
    constexpr GenericBitVector(GenericBitVector&& other) = default;
    
    // Move assignment OK
    GenericBitVector& operator=(GenericBitVector&& other) = default; 
    
private:
    
    // No copy
    GenericBitVector(const GenericBitVector& other) = delete;
    
    // No assignment
    GenericBitVector& operator=(const GenericBitVector& other) = delete;
    
    // Actual implementation
    // We need an encoder to actually make the bitvector.
    BitVectorEncoder* encoder;
    // And we need a place to put the vector when done.
    BitVector* bitvector;
    // And we need a persistent iterator
    BitVectorIterator* iterator;
    // And we need to remember how far along we are in our encoding.
    size_t size;
    
};

#endif
