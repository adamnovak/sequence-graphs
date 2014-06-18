// Test SmallSide objects.

#include "../SmallSide.hpp"

#include "SmallSideTests.hpp"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( SmallSideTests );

void SmallSideTests::setUp() {
}


void SmallSideTests::tearDown() {
}

/**
 * Test encoding and decoding SmallSides.
 */
void SmallSideTests::testEncodeDecode() {
    
    for(size_t id = 0; id < 1000; id++) {
        for(int face = 0; face < 2; face++) {
            // Construct a SmallSide
            SmallSide side(id, face);
            
            // Make sure we decode what we encoded
            CPPUNIT_ASSERT(side.getCoordinate() == id);
            CPPUNIT_ASSERT(side.getFace() == face);
        }
    }
}
