#ifndef GENERICBITVECTOR_HPP
#define GENERICBITVECTOR_HPP

#include <string>
#include <fstream>
#include <istream>
#include <ostream>
#include <utility>
#include "BitVector.hpp"

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
     * Add a 1 bit at the given index. Must be called in increasing order.
     * Cannot be run on a bitvector that has been loaded or finished.
     */
    void addBit(size_t index);
    
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
     * Get the size of the bitvector in bits. Must have been finished first.
     */
    size_t getSize() const;
    
    /**
     * Get the number of 1s occurring before the given index. Must be thread-
     * safe.
     */
    size_t rank(size_t index) const;
    
    /**
     * Get the number of 1s occurring before the given index. Must be thread-
     * safe. If at-least is set, includes a 1 at the given index.
     */
    size_t rank(size_t index, bool atLeast) const;
    
    /**
     * Returns true if the given index is set, or false otherwise. Must be
     * thread-safe.
     */
    bool isSet(size_t index) const;
    
    /**
     * Select the index at which the one with the given rank appears (i.e. the
     * first one would be 0). Must be thread-safe.
     */
    size_t select(size_t one) const;
    
    /**
     * OR two bitvectors together.
     */
    GenericBitVector* createUnion(const GenericBitVector& other) const;
    
    /**
     * Given an indes, return the index of the next 1 at or after that position,
     * paired with its rank.
     */
    std::pair<size_t, size_t> valueAfter(size_t index) const;
    
    /**
     * Given an indes, return the index of the last 1 at or before that
     * position, paired with its rank. Wraps around if no such value is found.
     */
    std::pair<size_t, size_t> valueBefore(size_t index) const;
    
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
    // We need an encoder to actually make the bitvector. We own this.
    BitVectorEncoder* encoder;
    // And we need a place to put the vector when done.
    BitVector* bitvector;
    // And we need to remember how far along we are in our encoding.
    size_t size;
    
};

#endif
