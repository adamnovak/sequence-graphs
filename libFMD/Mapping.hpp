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
    Mapping(TextPosition location, bool is_mapped=true, size_t context=0);
    
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
    inline TextPosition getLocation() {
        return location;
    }
    
    /**
     * Is this Mapping actually mapped?
     */
    inline bool isMapped() {
        return is_mapped;
    }
    
    /**
     * Return the amount of context used to map. This counts the base itself.
     */
    inline size_t getContext() {
        return context; 
    }
    
    // Holds (text, offset)
    TextPosition location;
    // Is the above actually filled in?
    bool is_mapped;
    // How much context was used for this?    
    size_t context;
};

/**
 * Provide pretty-printing for Mappings. See
 * <http://www.parashift.com/c++-faq/output-operator.html>
 */
std::ostream& operator<< (std::ostream& o, Mapping const& mapping);

#endif
