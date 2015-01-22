#include "MappingScheme.hpp"

MappingScheme::MappingScheme(const FMDIndex& index,
    const GenericBitVector& ranges, const GenericBitVector* mask): index(index),
    ranges(ranges), mask(mask), stats() {
    
    // Nothing to do!
}

StatTracker MappingScheme::getStats() const {
    return stats;
}
