 #include "DisambiguateFilter.hpp"
 #include <stdexcept>
 #include "Mapping.hpp"
 
 DisambiguateFilter::DisambiguateFilter(const FMDIndex& index): index(index) {
    // Nothing to do!
 }
 
 std::vector<Mapping> DisambiguateFilter::apply(
    std::vector<Mapping> leftMappings, std::vector<Mapping> rightMappings) {

    // This will hold all the disambiguated positions.
    std::vector<Mapping> toReturn;
    
    for(size_t i = 0; i < leftMappings.size(); i++) {
        // Disambiguate the pair of mappings, given that they both have the same
        // semantics.
        toReturn.push_back(index.disambiguate(leftMappings[i],
            rightMappings[i]));
    }
    
    // Return the vector of disambiguated mappings.
    return std::move(toReturn);
    
}
