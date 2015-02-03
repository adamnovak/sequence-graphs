#ifndef MATCHING_HPP
#define MATCHING_HPP

#include <iostream>
#include <cstdlib>
#include <functional>

#include "TextPosition.hpp"

/**
 * Keep track of an occurrence of a unique string which matches a query range to
 * a reference range. Since the query is implied, a Matching is defined by a
 * start, a length, and a location to which it matches.
 */
struct Matching {
    /**
     * Where does it start in the query?
     */
    size_t start;
    
    /**
     * How long is it?
     */
    size_t length;
    
    /**
     * Where is it in the reference?
     */
    TextPosition location;
    
    /**
     * Is this Matching less than that one of the same type on the same
     * query?
     */
    inline bool operator<(const Matching& other) const {
        // Order by start, then by length.
        return (start < other.start) ||
            (start == other.start && length < other.length);
    }
    
    /**
     * Determine if this Matching is equal to another.
     */
    inline bool operator==(const Matching& other) const {
        return start == other.start && length == other.length &&
            location == other.location; 
    }
    
    /**
     * Determine if this Matching is not equal to another.
     */
    inline bool operator!=(const Matching& other) const {
        return !(*this == other); 
    }
    
    /**
     * Make a new Matching.
     */
    inline Matching(size_t start, TextPosition location, size_t length): 
        start(start), location(location), length(length) {
    }
};

/**
 * Write this Matching to an output stream. Reports query coordinates.
 */
inline std::ostream& operator<<(std::ostream& out, const Matching& thing) {

    return out << thing.start << " - " << thing.start + thing.length;
}

namespace std {

    /**
     * Define a hash function for Matchings so they can be used as keys in an
     * std::unordered_map.
     */
    template <>
    struct hash<Matching> {

        /**
         * Hash a Matching to a size_t.
         */
        inline size_t operator()(const Matching& key) const {

            // We need a way to hash size_ts.
            std::hash<size_t> hasher;
            
            // Start with the start position
            size_t hash = hasher(key.start);
            
            // Mix in the length
            hash <<= 1;
            hash ^= hasher(key.length);
            // And the text
            hash <<= 1;
            hash ^= hasher(key.location.getText());
            // And the offset
            hash <<= 1;
            hash ^= hasher(key.location.getOffset());
            
            // And return the hash.
            return hash;
        }
    };

}
    
#endif
