#include "ZipMappingScheme.hpp"
#include "Log.hpp"
#include "util.hpp"

#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>

#include <boost/range/adaptor/reversed.hpp>

// Only template specializations can really be here.

template<>
bool ZipMappingScheme<FMDPosition>::canExtendThrough(FMDPosition context,
    const std::string& opposingQuery) const {
    
    Log::debug() << "Trying to extend " << context << " through " << 
        opposingQuery.size() << " opposing context" << std::endl;
    stats.add("extendThroughAttampts", 1);
        
    // We're going to retract it until it's no longer unique, then go
    // back and retract it one less.
    FMDPosition noLongerUnique = context;
    size_t nonUniqueLength = view.getIndex().retractRightOnly(noLongerUnique);
    
    while(noLongerUnique.isUnique(view)) {
        // It may take multiple retracts because of masks and stuff.
        nonUniqueLength = view.getIndex().retractRightOnly(noLongerUnique);
    }
    
    // Go back and retract one less (to a length one longer).
    FMDPosition barelyUnique = context;
    view.getIndex().retractRightOnly(barelyUnique, nonUniqueLength + 1);

    for(size_t i = opposingQuery.size() - 2; i != (size_t) -1 &&
        !barelyUnique.isEmpty(view); i--) {
        // Now go extend through with the opposing string, from right to left.
        // We skip the rightmost character (the one we are actually in the
        // process of mapping) because that's already searched, and we keep
        // going until we run out of results or we make it all the way through.
        
        Log::trace() << "Extending with " << opposingQuery[i] << std::endl;
        
        view.getIndex().extendLeftOnly(barelyUnique, opposingQuery[i]);
    }
    
    if(!barelyUnique.isEmpty(view)) {
        // We extended through!
        stats.add("extendThroughSuccesses", 1);
        return true;
    } else {
        // We didn't get any results upon extending through
        Log::debug() << "Extension failed" << std::endl;
        return false;
    }
}



