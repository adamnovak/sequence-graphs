#ifndef LRMAPPINGSCHEMETESTS_HPP
#define LRMAPPINGSCHEMETESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../LRMappingScheme.hpp"

/**
 * Mismatch Test for the FMDIndex.
 */
class LRMappingSchemeTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(LRMappingSchemeTests);
    CPPUNIT_TEST(testMap);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
    // Keep a pointer to a single index. That index is const, so we can't mess
    // it up between test cases.
    FMDIndex const* index;
    
    // Keep an LRMappingScheme around.
    LRMappingScheme* scheme;
    
    // And a ranges bitvector of all 1s
    GenericBitVector* ranges;
    
public:
    LRMappingSchemeTests();
    ~LRMappingSchemeTests();
    
    void setUp();
    void tearDown();

    void testMap();
    void testMapWithMask();
    
};

#endif
