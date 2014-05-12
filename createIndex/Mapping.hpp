#ifndef MAPPING_HPP
#define MAPPING_HPP

#include <iostream>
#include <utility>

#include "TextPosition.hpp"

/**
 * Represents a mapping between a base in a string and a (text, index) position
 * in the FMD-index. Contains the text and offset to which a character maps, and
 * a flag to say if it represents a real mapping or a result of "unmapped".
 */
struct Mapping
{
    // Holds (text, offset)
    TextPosition location;
    bool is_mapped;
    Mapping();
    Mapping(std::pair<int64_t, int64_t> location, bool is_mapped=true);
    /**
     * Provide equality comparison for testing.
     */
    bool operator==(const Mapping& other) const;
};

/**
 * Provide pretty-printing for Mappings. See
 * <http://www.parashift.com/c++-faq/output-operator.html>
 */
std::ostream& operator<< (std::ostream& o, Mapping const& mapping);

#endif
