#ifndef BWTTESTS_HPP
#define BWTTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

/**
 * Tests for the BWT.
 */
class BWTTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(BWTTests);
    CPPUNIT_TEST(testBWT);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;

public:
    void setUp();
    void tearDown();

    void testBWT();
};

#endif
