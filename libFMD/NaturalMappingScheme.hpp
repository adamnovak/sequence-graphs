#ifndef NATURALMAPPINGSCHEME_HPP
#define NATURALMAPPINGSCHEME_HPP

#include "MappingScheme.hpp"

/**
 * Represents a natural mapping scheme with only inexact credit support.
 */
class NaturalMappingScheme: public MappingScheme {

public:
    /**
     * Inherit the constructor.
     */
    using MappingScheme::MappingScheme;
    
    /**
     * Map the given query string according to the natural mapping scheme.
     * When a mapping is found, the callback function will be called with the
     * query base index, and the TextPosition to which it maps in the forward
     * direction.
     *
     */
    virtual void map(const std::string& query,
        std::function<void(size_t, TextPosition)> callback) const override;
        
    
    // Now come the scheme parameters and their default values.
    
    /**
     * How many mismatches can be tolerated when applying credit?
     */
    size_t z_max = 0;
    
    /**
     * Should credit be used to create more mappings?
     */
    bool credit = false;
    
    /**
     * What's the minimum context length to match on?
     */
    size_t minContext = 0;
    
};
   
#endif
