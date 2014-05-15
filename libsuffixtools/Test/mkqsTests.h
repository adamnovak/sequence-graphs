#ifndef MKQSTESTS_HPP
#define MKQSTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

/**
 * Tests for multi-key quicksort.
 */
class MKQSTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(MKQSTests);
    CPPUNIT_TEST(testSort);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;

private:
    ReadTable* readTable;
    SuffixArray* suffixArray;

public:
    void setUp();
    void tearDown();

    void testSort();
};

#endif
