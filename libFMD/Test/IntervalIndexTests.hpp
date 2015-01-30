#ifndef INTERVALINDEXTESTS_HPP
#define INTERVALINDEXTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

#include <string>

/**
 * Tests for SmallSide.
 */
class IntervalIndexTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(IntervalIndexTests);
    CPPUNIT_TEST(testCreate);
    CPPUNIT_TEST(testLookupStartingBefore);
    CPPUNIT_TEST(testLookupStartingAfter);
    CPPUNIT_TEST(testLookupEndingAfter);
    CPPUNIT_TEST(testLookupEndingBefore);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void setUp();
    void tearDown();

    void testCreate();
    void testLookupStartingBefore();
    void testLookupStartingAfter();
    void testLookupEndingAfter();
    void testLookupEndingBefore();
};

#endif
