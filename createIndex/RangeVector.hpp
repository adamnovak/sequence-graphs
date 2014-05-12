#ifndef RANGEVECTOR_HPP
#define RANGEVECTOR_HPP

/**
 * Defines a RangeVector type that is a bit vector we use for ranges, with
 * efficient rank and select.
 */

#ifdef USE_NIBBLE_VECTORS
    // Use Nibble Vectors to encode our range endpoint bitmaps
    #include "rlcsa/bits/nibblevector.h"
    typedef CSA::NibbleVector RangeVector;
    typedef CSA::NibbleEncoder RangeEncoder;
#else
    // Use RLEVectors to encode our range endpoint bitmaps
    #include "rlcsa/bits/rlevector.h"
    typedef CSA::RLEVector RangeVector;
    typedef CSA::RLEEncoder RangeEncoder;
#endif
 

#endif
