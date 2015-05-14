#ifndef MAPPINGSCHEME_HPP
#define MAPPINGSCHEME_HPP

#include "GenericBitVector.hpp"
#include "FMDIndexView.hpp"
#include "TextPosition.hpp"
#include "StatTracker.hpp"

#include <string>
#include <functional>
#include <map>

/**
 * Represents a mapping algorithm and its associated parameters. Subclasses
 * actually implement it.
 */
class MappingScheme {
public:

    /**
     * Make a new MappingScheme against the given view of an FMDIndex. Takes the
     * view for itself.
     *
     * This constructor ought to be inherited by all subclasses with "using".
     * All other parameters of subclasses ought to be optional, with sensible
     * default values, and settable by public field access.
     */
    MappingScheme(FMDIndexView&& index);

    // Default copy/move constructors OK    
    MappingScheme(const MappingScheme& other) = default;
    MappingScheme(MappingScheme&& other) = default;
    
    // Assignment not OK, since we use references. But since we need
    // polymorphism nobody will ever be able to use us by value anyway...
    MappingScheme& operator=(const MappingScheme& other) = delete;
    MappingScheme& operator=(MappingScheme&& other) = delete;

    // Defaut destructor also OK.
    
    /**
     * Map the given query string according to the mapping algorithm. When a
     * mapping is found, the callback function will be called with the query
     * base index, and the TextPosition to which it maps in the forward
     * direction.
     *
     * Must be defined by all implementations.
     *
     * Must be thread-safe.
     *
     * Should update the "mapped" and "unmapped" stats in the MappingScheme's
     * StatTracker, as well as any other applicable stats.
     */
    virtual void map(const std::string& query,
        std::function<void(size_t, TextPosition)> callback) const = 0;
    
    /**
     * Get a snapshot of the stats for this mapping scheme.
     */
    StatTracker getStats() const;
    
protected:
    // These configuration parameters are ones that every MappingScheme that
    // uses an FMDIndex (all the ones we care about) will need.
    
    /**
     * The FMDIndexView against which we are mapping. The FMDIndex which we map
     * against is accessible through this.
     */
    const FMDIndexView view;
    
    /**
     * We keep a StatTracker around for tracking stats. It is mutable so that
     * the map method can update it, and it is guaranteed to be thread-safe.
     */
    mutable StatTracker stats;
};
 

#endif
