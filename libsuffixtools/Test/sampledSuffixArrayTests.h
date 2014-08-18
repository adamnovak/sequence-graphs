#ifndef SAMPLEDSUFFIXARRAYTESTS_HPP
#define SAMPLEDSUFFIXARRAYTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

/**
 * Tests for Sampled Suffix Array generartion.
 */
class SampledSuffixArrayTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(SampledSuffixArrayTests);
    CPPUNIT_TEST(testConstruction);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;

private:
    // Test fixtures
    ReadTable* readTable;
    ReadInfoTable* infoTable;
    SuffixArray* suffixArray;

public:
    void setUp();
    void tearDown();

    void testConstruction();
};

#endif
