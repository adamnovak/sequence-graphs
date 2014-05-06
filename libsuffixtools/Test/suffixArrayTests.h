#ifndef SUFFIXARRAYTESTS_HPP
#define SUFFIXARRAYTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

/**
 * Tests for the SuffixArray.
 */
class SuffixArrayTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(SuffixArrayTests);
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
