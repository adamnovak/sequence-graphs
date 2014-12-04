#ifndef MAPPINGSCHEME_HPP
#define MAPPINGSCHEME_HPP

#include "GenericBitVector.hpp"
#include "FMDIndex.hpp"
#include "TextPosition.hpp"

#include <string>
#include <functional>

/**
 * Represents a mapping algorithm and its associated parameters. Subclasses
 * actually implement it.
 */
class MappingScheme {
public:

    /**
     * Make a new MappingScheme against the given FMDIndex, mapping to the given
     * ranges, and using the given mask of mappable BWT positions.
     *
     * This constructor ought to be inherited by all subclasses with "using".
     * All other parameters of subclasses ought to be optional, with sensible
     * default values, and settable by public field access.
     */
    MappingScheme(const FMDIndex& index, const GenericBitVector& ranges,
        const GenericBitVector* mask = NULL);

    
    /**
     * Map the given query string according to the mapping algorithm. When a
     * mapping is found, the callback function will be called with the query
     * base index, and the TextPosition to which it maps in the forward
     * direction.
     *
     * Must be defined by all implementations.
     */
    virtual void map(const std::string& query,
        std::function<void(size_t, TextPosition)> callback) = 0;
    
protected:
    // These configuration parameters are ones that every MappingScheme that
    // uses an FMDIndex (all the ones we care about) will need.
    
    /**
     * The FMDIndex against which we are mapping.
     */
    const FMDIndex& index;
    
    /**
     * The bit vector defining ranges to map to on the index.
     */
    const GenericBitVector& ranges;
    
    /**
     * The mask with which we are mapping.
     */
    const GenericBitVector* mask;
}
 

#endif
