// Test GenericBitVector objects.

#include "../GenericBitVector.hpp"

#include "GenericBitVectorTests.hpp"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( GenericBitVectorTests );

void GenericBitVectorTests::setUp() {
}


void GenericBitVectorTests::tearDown() {
}


/**
 * Make sure isSet works.
 */
void GenericBitVectorTests::testIsSet() {

    // Make a pair of encoders.
    GenericBitVector v;
    BitVectorEncoder v2(32);
    
    // Populate them
    v.addBit(1);
    v2.addBit(1);
    
    v.addBit(3);
    v2.addBit(3);
    
    v.addBit(5);
    v2.addBit(5);
    
    // Finish up the new one
    v.finish(10);
    
    // And the old one
    v2.flush();
    BitVector v2vector(v2, 10);
    BitVectorIterator v2iterator(v2vector);
    
    for(size_t i = 0; i < 20; i++) {
        // Check a bunch of positions.
        CPPUNIT_ASSERT(v.isSet(i) == v2iterator.isSet(i));
    }
    
    
}


void GenericBitVectorTests::testRank() {
}
void GenericBitVectorTests::testSelect() {
}
void GenericBitVectorTests::testValueBefore() {
}
void GenericBitVectorTests::testValueAfter() {
}

