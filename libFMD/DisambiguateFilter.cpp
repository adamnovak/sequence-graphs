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
                
                // Flip the right one.
                TextPosition rightPosition = rightMappings[i].getLocation();
                rightPosition.flip(index.getContigLength(index.getContigNumber(
                    rightPosition)));
                
                if(leftMappings[i].getLocation() == rightPosition) {
                    // They agree. Take either.
                    toReturn.push_back(Mapping(rightPosition));
                } else {
                    // They disagree. Say unmapped.
                    toReturn.push_back(Mapping());
                }
            } else {
                // Only left is mapped. Take it.
                toReturn.push_back(leftMappings[i].getLocation());
            }
        } else if(rightMappings[i].isMapped()) {
            // Flip the right mapping and take that
            TextPosition rightPosition = rightMappings[i].getLocation();
            rightPosition.flip(index.getContigLength(index.getContigNumber(
                rightPosition)));
                
            toReturn.push_back(Mapping(rightPosition));
            
        } else {
            // No mapping, say unmapped.
            toReturn.push_back(Mapping());
        }
    }
    
    return std::move(toReturn);
    
}
