#ifndef FMDINDEXBUILDERTESTS_HPP
#define FMDINDEXBUILDERTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "FMDIndex.hpp"

/**
 * Tests for the FMDIndex.
 */
class FMDIndexTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(FMDIndexTests);
    CPPUNIT_TEST(testDump);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
public:
    void setUp();
    void tearDown();

    void testDump();
};

#endif
