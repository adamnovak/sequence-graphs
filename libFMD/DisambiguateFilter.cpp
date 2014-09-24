 #include "DisambiguateFilter.hpp"
 
 DisambiguateFilter::DisambiguateFilter(FMDIndex& index): index(index) {
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
                    toReturn.append(Mapping(rightPosition));
                } else {
                    // They disagree. Say unmapped.
                    toReturn.append(Mapping());
                }
            } else {
                // Only left is mapped. Take it.
                toReturn.append(leftMappings[i].getLocation());
            }
        } else if(rightMappings[i].isMapped()) {
            // Flip the right mapping and take that
            TextPosition rightPosition = rightMappings[i].getLocation();
            rightPosition.flip(index.getContigLength(index.getContigNumber(
                rightPosition)));
                
            toReturn.append(Mapping(rightPosition))
            
        } else {
            // No mapping, say unmapped.
            toReturn.append(Mapping());
        }
    }
    
    return std::move(toReturn);
    
}
