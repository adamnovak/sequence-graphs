#ifndef NATURALMAPPINGSCHEME_HPP
#define NATURALMAPPINGSCHEME_HPP

#include "MappingScheme.hpp"
#include "Mapping.hpp"

/**
 * Mapping scheme implementing Benedict's "natural" mapping scheme. If all
 * the unique-in-the-reference strings overlapping a query base agree about
 * where the abse should be, and you can't get from the base off the end of
 * the query with multiple results, map the base.
 */
class NaturalMappingScheme: public MappingScheme {

public:
    /**
     * Inherit the constructor.
     */
    using MappingScheme::MappingScheme;
    
    /**
     * Map the given query string according to the natural mapping scheme with
     * optional inexact credit. When a mapping is found, the callback function
     * will be called with the query base index, and the TextPosition to which
     * it maps in the forward direction.
     *
     */
    virtual void map(const std::string& query,
        std::function<void(size_t, TextPosition)> callback) const override;
        
    
    // Now come the scheme parameters and their default values.
    
    /**
     * What's the minimum factor by which the maximum unique string length must
     * exceed the minimum?
     */
    double multContext = 0.0;
    
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
    
protected:
    /**
     * Map the given query string, producing a vector of Mappings. Does not
     * include credit yet.
     */
    std::vector<Mapping> naturalMap(const std::string& query) const;
    
    /**
     * Keep track of an occurrence of a unique string which matches a query
     * position to a reference position.
     */
    struct Matching {
        /**
         * Where does it start in the query?
         */
        size_t start;
        /**
         * Where is it in the reference?
         */
        TextPosition location;
        /**
         * How long is it?
         */
        size_t length;
        
        /**
         * Make a new Matching.
         */
        inline Matching(size_t start, TextPosition location, size_t length): 
            start(start), location(location), length(length) {
        };
    };
    
    /**
     * Find all of the maximal unique matchings between query string characters
     * and the reference.
     */
    std::vector<Matching> findMaxMatchings(const std::string& query) const;
        
    /**
     * Find all of the minimal unique matchings between query string characters
     * and the reference.
     */
    std::vector<Matching> findMinMatchings(const std::string& query) const;
    
};
   
#endif
