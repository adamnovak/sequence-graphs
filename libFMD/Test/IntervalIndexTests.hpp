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
    CPPUNIT_TEST(testLookupBefore);
    CPPUNIT_TEST(testLookupAfter);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void setUp();
    void tearDown();

    void testCreate();
    void testLookupBefore();
    void testLookupAfter();
};

#endif
