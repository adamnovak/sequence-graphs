#ifndef FMDINDEXMISMATCHTESTS_HPP
#define FMDINDEXMISMATCHTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../FMDIndex.hpp"

/**
 * Mismatch Test for the FMDIndex.
 */
class FMDIndexMismatchTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(FMDIndexMismatchTests);
    CPPUNIT_TEST(testMismatch);
    CPPUNIT_TEST(testMapOnMismatch);
    CPPUNIT_TEST(testMismatchCount);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
    // Keep a pointer to a single index. That index is const, so we can't mess
    // it up between test cases.
    FMDIndex const* index;
    
public:
    FMDIndexMismatchTests();
    ~FMDIndexMismatchTests();
    
    void setUp();
    void tearDown();

    void testMismatch();
    void testMapOnMismatch();
    void testMismatchCount();
    
};

#endif
