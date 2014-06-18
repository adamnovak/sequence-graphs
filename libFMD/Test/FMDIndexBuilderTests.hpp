#ifndef FMDINDEXBUILDERTESTS_HPP
#define FMDINDEXBUILDERTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../FMDIndexBuilder.hpp"

/**
 * Tests for the FMDIndexBuilder.
 */
class FMDIndexBuilderTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(FMDIndexBuilderTests);
    CPPUNIT_TEST(testBuild);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
public:
    void setUp();
    void tearDown();

    void testBuild();
    
};

#endif
