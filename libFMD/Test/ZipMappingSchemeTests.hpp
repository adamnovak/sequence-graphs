#ifndef ZIPMAPPINGSCHEMETESTS_HPP
#define ZIPMAPPINGSCHEMETESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../ZipMappingScheme.hpp"

/**
 * Test for the zip mapping scheme, using sequences we can merge.
 */
class ZipMappingSchemeTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(ZipMappingSchemeTests);
    CPPUNIT_TEST(testMap);
    CPPUNIT_TEST(testMapWithMask);
    CPPUNIT_TEST(testMapWithMaskAndRanges);
    CPPUNIT_TEST(testMapWithGroups);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
    // Keep a pointer to a single index. That index is const, so we can't mess
    // it up between test cases.
    FMDIndex const* index;
    
    // Keep a ZipMappingScheme around. We want to try both template
    // specializations, so we use the base class type for the pointer.
    MappingScheme* scheme;
    
    // And a ranges bitvector merging some positions
    GenericBitVector* ranges;
    
public:
    ZipMappingSchemeTests();
    ~ZipMappingSchemeTests();
    
    void setUp();
    void tearDown();

    void testMap();
    void testMapWithMask();
    void testMapWithMaskAndRanges();
    void testMapWithGroups();
};

#endif
