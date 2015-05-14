#include "MappingScheme.hpp"

MappingScheme::MappingScheme(FMDIndexView&& view): view(std::move(view)),
    stats() {

    // Nothing to do!
}

StatTracker MappingScheme::getStats() const {
    return stats;
}
