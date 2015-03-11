#include "MappingScheme.hpp"

MappingScheme::MappingScheme(const FMDIndexView view): view(std::move(view)),
    stats() {

    // Nothing to do!
}

StatTracker MappingScheme::getStats() const {
    return stats;
}
