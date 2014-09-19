#ifndef ADJACENCYCOMPONRNTUTIL_HPP
#define ADJACENCYCOMPONENTUTIL_HPP

#include <vector>
#include <cstdint> 
#include <map>
#include <functional>

// Grab pinchesAndCacti dependency.
#include <stPinchGraphs.h>

/**
 * adjacencyComponentUtil.hpp: utility functions for working with adjacency
 * components from a pinch graph.
 */
 
 /**
 * Take a pinch thread set and get the adjacency components in a C++ idiomatic
 * data structure.
 *
 * Returns a vector of adjacency components, which are vectors of pinch ends.
 */
std::vector<std::vector<stPinchEnd>> 
getAdjacencyComponents(
    stPinchThreadSet* threadSet
);

/**
 * Take a vector of adjacency components and get the spectrum of adjacency
 * component sizes. Size 2 components are things like SNPs and indels, while
 * larger components are probably more complex structures.
 */
std::map<size_t, size_t>
getAdjacencyComponentSpectrum(
    std::vector<std::vector<stPinchEnd>> components
);

/**
 * Get the components of a certain size from a vector of adjacency components.
 * Copies all those components.
 */
std::vector<std::vector<stPinchEnd>>
filterComponentsBySize(
    std::vector<std::vector<stPinchEnd>> components,
    size_t size,
    std::function<bool(size_t, size_t)> = [](size_t a, size_t b) { 
        return a == b;
    }
);

/**
 * Given a pair of pinch ends, get all the segments on all paths from the first
 * to the second. There must be no paths leaving the first on adjacency edges
 * that do not enter the second on an adjacency edge.
 */
std::vector<std::vector<stPinchSegment*>>
getAllPaths(
    stPinchEnd* start,
    stPinchEnd* end
);

/**
 * Take a vector of adjacency components all of size 2, and work out the indel
 * length difference for all of them that are simple indels. A straight
 * substitution is a length-0 indel.
 */
std::vector<int64_t>
getIndelLengths(
    std::vector<std::vector<stPinchEnd>> components
);

/**
 * Take a vector of adjacency components all of size 4, and count all the tandem
 * duplications. We call it a tandem duplication if two of the ends belong to
 * the same block and are connected.
 */
size_t
countTandemDuplications(
    std::vector<std::vector<stPinchEnd>> components
);

/**
 * Save the given list of adjacency components to the given filename.
 */
void
writeAdjacencyComponents(
    std::vector<std::vector<stPinchEnd>> components,
    const std::string& filename
);



#endif
