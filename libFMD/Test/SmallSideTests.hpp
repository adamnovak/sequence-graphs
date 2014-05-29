#ifndef SMALLSIDETESTS_HPP
#define SMALLSIDETESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

/**
 * Tests for SmallSide.
 */
class SmallSideTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(SmallSideTests);
    CPPUNIT_TEST(testEncodeDecode);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void setUp();
    void tearDown();

    void testEncodeDecode();
};

#endif
