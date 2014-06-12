#ifndef RANGEVECTOR_HPP
#define RANGEVECTOR_HPP

/**
 * Defines a BitVector type that is a bit vector we use for ranges and genome
 * masking, with efficient rank and select.
 */

// Use Nibble Vectors, borrowed from RLCSA, to encode our range endpoint bitmaps
#include "CSA/NibbleVector.hpp"
typedef CSA::NibbleVector BitVector;
typedef CSA::NibbleVector::Iterator BitVectorIterator;
typedef CSA::NibbleEncoder BitVectorEncoder;

 

#endif
