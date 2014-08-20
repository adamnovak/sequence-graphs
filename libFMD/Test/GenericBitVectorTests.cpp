// Test GenericBitVector objects.

#include "../GenericBitVector.hpp"

#include "GenericBitVectorTests.hpp"
#include "../Log.hpp"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( GenericBitVectorTests );

void GenericBitVectorTests::setUp() {
}


void GenericBitVectorTests::tearDown() {
}

/**
 * Get a pair of test fixtures.
 */
std::pair<GenericBitVector*, BitVector*> 
    GenericBitVectorTests::makeTestData() {

    // Make a pair of encoders.
    GenericBitVector* v = new GenericBitVector;
    BitVectorEncoder v2(32);
    
    // Populate them
    v->addBit(2);
    v2.addBit(2);
    
    v->addBit(4);
    v2.addBit(4);
    
    v->addBit(5);
    v2.addBit(5);
    
    v->addBit(6);
    v2.addBit(6);
    
    // Finish up the new one
    v->finish(10);
    
    // And the old one
    v2.flush();
    BitVector* v2vector = new BitVector(v2, 10);
    

    return std::make_pair(v, v2vector);
}

/**
 * Make sure isSet works.
 */
void GenericBitVectorTests::testIsSet() {

    auto pair = makeTestData();
    GenericBitVector* v = pair.first;
    BitVectorIterator* v2 = new BitVectorIterator(*pair.second);
    
    for(size_t i = 0; i < 20; i++) {
        // Check a bunch of positions.
        CPPUNIT_ASSERT(v->isSet(i) == v2->isSet(i));
    }
    
    delete pair.first;
    delete pair.second;
    delete v2;
    
}


void GenericBitVectorTests::testRank() {

    auto pair = makeTestData();
    GenericBitVector* v = pair.first;
    BitVectorIterator* v2 = new BitVectorIterator(*pair.second);
    
    for(size_t i = 0; i < 20; i++) {
        // Check a bunch of positions.
        CPPUNIT_ASSERT(v->rank(i) == v2->rank(i));
        CPPUNIT_ASSERT(v->rank(i, true) == v2->rank(i, true));
        CPPUNIT_ASSERT(v->rank(i, false) == v2->rank(i, false));
    }
    
    delete pair.first;
    delete pair.second;
    delete v2;

}
void GenericBitVectorTests::testSelect() {

    auto pair = makeTestData();
    GenericBitVector* v = pair.first;
    BitVectorIterator* v2 = new BitVectorIterator(*pair.second);
    
    for(size_t i = 0; i < 6; i++) {
        // Check a bunch of positions.
        CPPUNIT_ASSERT(v->select(i) == v2->select(i));
    }
    
    delete pair.first;
    delete pair.second;
    delete v2;

}
void GenericBitVectorTests::testValueBefore() {
    
    auto pair = makeTestData();
    GenericBitVector* v = pair.first;
    BitVectorIterator* v2 = new BitVectorIterator(*pair.second);
    
    for(size_t i = 0; i < 20; i++) {
        // Check a bunch of positions.
        auto ourResult = v->valueBefore(i);
        auto theirResult = v2->valueBefore(i);
        Log::debug() << "Position " << i << ": " << ourResult.first << ", " <<
            ourResult.second << " vs. " << theirResult.first << ", " << 
            theirResult.second << std::endl;
        CPPUNIT_ASSERT(ourResult == theirResult);
    }
    
    delete pair.first;
    delete pair.second;
    delete v2;
}
void GenericBitVectorTests::testValueAfter() {

    auto pair = makeTestData();
    GenericBitVector* v = pair.first;
    BitVectorIterator* v2 = new BitVectorIterator(*pair.second);
    
    for(size_t i = 0; i < 20; i++) {
        // Check a bunch of positions.
        auto ourResult = v->valueAfter(i);
        auto theirResult = v2->valueAfter(i);
        Log::debug() << "Position " << i << ": " << ourResult.first << ", " <<
            ourResult.second << " vs. " << theirResult.first << ", " << 
            theirResult.second << std::endl;
        CPPUNIT_ASSERT(ourResult == theirResult);
    }
    
    delete pair.first;
    delete pair.second;
    delete v2;
    
}

