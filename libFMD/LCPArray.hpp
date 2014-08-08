#ifndef LCPARRAY_HPP
#define LCPARRAY_HPP

#include <vector>
#include <string>

// Depend on the libsuffixtools stuff.
#include <ReadTable.h>
#include <SuffixArray.h>

#include "Log.hpp"

/**
 * Defines an array suitable for holding Longest Common Prefix information
 * between successive suffixes in a suffix array or FMD-index. Allows efficient
 * queries for the position of the next or previous smaller value from any
 * position.
 *
 * This particular implementation stores a full index, and so will use an
 * inordinately large amount of space.
 */
class LCPArray {

public:
    /**
     * Create a new LCPArray by building it from the given suffix array, using
     * the given collection of strings to get the actual characters.
     *
     * Neither the suffix array nor the read table are needed after the
     * constructor returns.
     */
    LCPArray(const SuffixArray& suffixArray, const ReadTable& strings);
    
    /**
     * Load an LCPArray from the given file. Uses platform-dependent byte
     * order and size_t size.
     */
    LCPArray(const std::string& filename);
    
    /**
     * Save an LCPArray to the given file. Uses platform-dependent byte order
     * and size_t size.
     */
    void save(const std::string& filename) const;
    
    /**
     * Get the longest common prefix value at a given index.
     */
    inline size_t operator[](size_t index) const {
        return values[index];
    }
    
    /**
     * Get the index of the previous smaller value before the given index.
     */
    inline size_t getPSV(size_t index) const {
        return psvs[index];
    }
    
    /**
     * Get the index of the next smaller value after the given index.
     */
    inline size_t getNSV(size_t index) const {
        return nsvs[index];
    }

protected:
    // TODO: Should these helpers go on the increasingly inaccurately named
    // ReadTable?

    /**
     * Get the character at the given offset into the given suffix in the given
     * collection of strings.
     */
    static inline char getFromSuffix(const SAElem& suffix, size_t offset,
        const ReadTable& strings) {
    
        // Split out the text and offset from the suffix, and offset the offset.
        return strings.getChar(suffix.getID(), suffix.getPos() + offset);    
        
    }
        
    /**
     * Get the length of a suffix form the collection of strings.
     */
    static inline size_t getSuffixLength(const SAElem& suffix,
        const ReadTable& strings) {
        
        // Get the length of the "read", minus the offset.
        return strings.getReadLength(suffix.getID()) - suffix.getPos();
        
    }
    
    /**
     * Increment the given SAElem to the next valid position in the hypothetical
     * single string.
     */
    static inline void increment(SAElem& target, const ReadTable& strings) {
        Log::trace() << "Incrementing" << std::endl;
        
        // How much of this suffix is left?
        size_t suffixLength = getSuffixLength(target, strings);
        
        if(suffixLength > 0) {
            // There's room to advance in the string
            target.setPos(target.getPos() + 1);
        } else {
            // We need to advance to the start of the next string
            target.setID(target.getID() + 1);
            target.setPos(0);
        }
    }
    
    /**
     * Increment the given SAElem by a given amount in the hypothetical string.
     */
    static inline void incrementBy(SAElem& target, size_t offset,
        const ReadTable& strings) {
        
        Log::trace() << "Incrementing by " << offset << std::endl;
        
        while(offset > 0) {
            
            Log::trace() << "Offset: " << offset << std::endl;
        
            // How much of this suffix is left?
            size_t suffixLength = getSuffixLength(target, strings);
            
            if(suffixLength > offset) {
                // There's room to advance in the string. Just do that.
                target.setPos(target.getPos() + offset);
                offset = 0;
            } else {
                // Go to the next string and see how much room that buys us.
                offset -= suffixLength;
                target.setID(target.getID() + 1);
                target.setPos(0);
            }
        }
        
    }
    
    /**
     * Returns true if the given suffix has not fallen off the end of the
     * string, false otherwise.
     */
    static inline bool inRange(const SAElem& suffix, const ReadTable& strings) {
        return suffix.getID() < strings.getCount() && 
            getSuffixLength(suffix, strings) > 0;
    }
    
    // Store all the LCP array entries.
    std::vector<size_t> values;
    
    // Store the index of the previous smaller value for each position.
    std::vector<size_t> psvs;
    
    // Store the index of the next smaller value for each position.
    std::vector<size_t> nsvs;
    
};

#endif
