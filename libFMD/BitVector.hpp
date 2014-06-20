#ifndef RANGEVECTOR_HPP
#define RANGEVECTOR_HPP

/**
 * Defines a BitVector type that is a bit vector we use for ranges and genome
 * masking, with efficient rank and select.
 */

// Use Nibble Vectors, borrowed from RLCSA, and renamed to BitVectors, to encode
// our range endpoint bitmaps
#include "CSA/BitVector.hpp"

// Previously we just typedef'd NibbleVector thigns to BitVector things, but
// that confused SWIG.

using CSA::BitVector;
using CSA::BitVectorIterator;
using CSA::BitVectorEncoder;
 

#endif
