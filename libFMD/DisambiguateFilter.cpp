 #include "DisambiguateFilter.hpp"
 
 DisambiguateFilter::DisambiguateFilter(const FMDIndex& index): index(index) {
    // Nothing to do!
 }
 
 std::vector<Mapping> DisambiguateFilter::apply(
    std::vector<Mapping> leftMappings, std::vector<Mapping> rightMappings) {

    // This will hold all the disambiguated positions.
    std::vector<Mapping> toReturn;
    
    for(size_t i = 0; i < leftMappings.size(); i++) {
        if(leftMappings[i].isMapped()) {
            if(rightMappings[i].isMapped()) {
                // Both are mapped.
                
                // Flip the left one.
                TextPosition leftPosition = leftMappings[i].getLocation();
                leftPosition.flip(index.getContigLength(index.getContigNumber(
                    leftPosition)));
                
                if(rightMappings[i].getLocation() == leftPosition) {
                    // They agree. Take either.
                    toReturn.push_back(Mapping(leftPosition));
                } else {
                    // They disagree. Say unmapped.
                    toReturn.push_back(Mapping());
                }
            } else {
                // Only left is mapped. Flip it and take it.
                TextPosition leftPosition = leftMappings[i].getLocation();
                leftPosition.flip(index.getContigLength(index.getContigNumber(
                    leftPosition)));
                toReturn.push_back(leftPosition);
            }
        } else if(rightMappings[i].isMapped()) {
            // Take the right mapping.
            toReturn.push_back(rightMappings[i]);
            
        } else {
            // No mapping, say unmapped.
            toReturn.push_back(Mapping());
        }
    }
    
    return std::move(toReturn);
    
}
