#ifndef GENERICBITVECTORTESTS_HPP
#define GENERICBITVECTORTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <utility>

#include "../GenericBitVector.hpp"
#include "../BitVector.hpp"

/**
 * Tests for GenericBitVector against BitVector.
 */
class GenericBitVectorTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(GenericBitVectorTests);
    CPPUNIT_TEST(testIsSet);
    CPPUNIT_TEST(testRank);
    CPPUNIT_TEST(testSelect);
    CPPUNIT_TEST(testValueBefore);
    CPPUNIT_TEST(testValueAfter);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void setUp();
    void tearDown();

    void testIsSet();
    void testRank();
    void testSelect();
    void testValueBefore();
    void testValueAfter();
    
    std::pair<GenericBitVector*, BitVector*> makeTestData();
};

#endif
