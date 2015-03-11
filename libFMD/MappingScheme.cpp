#include "MappingScheme.hpp"

MappingScheme::MappingScheme(const FMDIndexView& view): view(view), stats() {
    // Nothing to do!
}

StatTracker MappingScheme::getStats() const {
    return stats;
}
