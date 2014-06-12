#ifndef RANGEVECTOR_HPP
#define RANGEVECTOR_HPP

/**
 * Defines a RangeVector type that is a bit vector we use for ranges, with
 * efficient rank and select.
 */

// Use Nibble Vectors, borrowed from RLCSA, to encode our range endpoint bitmaps
#include "CSA/NibbleVector.hpp"
typedef CSA::NibbleVector RangeVector;
typedef CSA::NibbleVector::Iterator RangeVectorIterator;
typedef CSA::NibbleEncoder RangeEncoder;

 

#endif
