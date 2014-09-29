#ifndef MAPPING_HPP
#define MAPPING_HPP

#include <iostream>
#include <utility>

#include "TextPosition.hpp"

/**
 * Represents a mapping between a base in a string and a (text, index) position
 * in the FMD-index. Contains the text and offset to which a character maps, and
 * a flag to say if it represents a real mapping or a result of "unmapped". Also
 * contains how much (maximal) context was used to map this position. (Generally
 * this should be the maximum possiblecontext for credit to work correctly.)
 */
struct Mapping
{
public:
    Mapping();
    Mapping(TextPosition location);
    Mapping(TextPosition location, size_t leftContext, size_t rightContext);
    
    /**
     * Provide equality comparison for testing.
     */
    bool operator==(const Mapping& other) const;
    
    /**
     * Also provide inequality.
     */
    bool operator!=(const Mapping& other) const;
    
    /**
     * What text and offset is this mapping to?
     */
    inline TextPosition getLocation() const {
        return location;
    }
    
    /**
     * Is this Mapping actually mapped?
     */
    inline bool isMapped() const {
        return is_mapped;
    }
    
    /**
     * Return the amount of context used to map on the left. This counts the
     * base itself. Returns 0 if not mapped using this side. 0 for
     * unmapped Mappings.
     */
    inline size_t getLeftContext() const {
        return leftContext; 
    }
    
    /**
     * Return the amount of context used to map on the right. This counts the
     * base itself. Returns 0 if not mapped using this side. 0 for
     * unmapped Mappings.
     */
    inline size_t getRightContext() const {
        return rightContext; 
    }
    
    
    /**
     * Return the amount of total context used to map this base.
     */
    inline size_t getContext() const {
        // Don't double-count the central base if both contexts are nonzero.
        return leftContext + (rightContext > 0) ? rightContext - 1 : 0;
    }
    
    /**
     * Flip this mapping and return a new mapping for the other strand, with
     * contexts swapped.
     */
    inline Mapping flip(size_t contigLength) const {
        if(is_mapped) {
            // If we're mapped, flip the location and exchange the contexts.
            TextPosition newLocation = location;
            newLocation.flip(contigLength);
            return Mapping(newLocation, rightContext, leftContext);
        } else {
            // If we are unmapped, just continue being exactly the same.
            return *this;
        }
    }
    
    // Holds (text, offset)
    TextPosition location;
    // Is the above actually filled in?
    bool is_mapped;
    // How much left context was used for this?    
    size_t leftContext;
    // And how much right context?
    size_t rightContext;
};

/**
 * Provide pretty-printing for Mappings. See
 * <http://www.parashift.com/c++-faq/output-operator.html>
 */
std::ostream& operator<< (std::ostream& o, Mapping const& mapping);

#endif
