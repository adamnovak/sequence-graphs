#ifndef DISAMBIGUATEFILTERTESTS_HPP
#define DISAMBIGUATEFILTERTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../DisambiguateFilter.hpp"

/**
 * Tests for the DiasmbiguateFilter.
 */
class DisambiguateFilterTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(DisambiguateFilterTests);
    CPPUNIT_TEST(testApply);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
    // Keep a pointer to a single index. That index is const, so we can't mess
    // it up between test cases.
    FMDIndex const* index;
    
public:
    DisambiguateFilterTests();
    ~DisambiguateFilterTests();
    
    void setUp();
    void tearDown();

    void testApply();
    
};

#endif
