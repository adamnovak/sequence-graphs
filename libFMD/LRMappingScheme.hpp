#ifndef LRMAPPINGSCHEME_HPP
#define LRMAPPINGSCHEME_HPP

#include "MappingScheme.hpp"

/**
 * Represents a left-right mapping scheme supporting inexact matching and a
 * variety of context limits.
 */
class LRMappingScheme: public MappingScheme {

public:
    /**
     * Inherit the constructor.
     */
    using MappingScheme::MappingScheme;
    
    /**
     * Map the given query string according to the left-right mapping scheme.
     * When a mapping is found, the callback function will be called with the
     * query base index, and the TextPosition to which it maps in the forward
     * direction.
     *
     */
    virtual void map(const std::string& query,
        std::function<void(size_t, TextPosition)> callback) const override;
        
    
    // Now come the scheme parameters and their default values.
    
    /**
     * How many mismatches can be tolerated when mapping?
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
    
    /**
     * What's the minimum additional context to use after becoming unique?
     */
    size_t addContext = 0;
    
    /**
     * What's the minimum factor by which the maximum context length must exceed
     * the minimum?
     */
    double multContext = 0.0;
    
};
   
#endif
