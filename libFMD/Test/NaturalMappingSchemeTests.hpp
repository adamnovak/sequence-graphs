#ifndef NATURALMAPPINGSCHEMETESTS_HPP
#define NATURALMAPPINGSCHEMETESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../NaturalMappingScheme.hpp"

/**
 * Test for the natural mapping scheme.
 */
class NaturalMappingSchemeTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(NaturalMappingSchemeTests);
    CPPUNIT_TEST(testMap);
    CPPUNIT_TEST(testMapWithMask);
    CPPUNIT_TEST(testSkipMismatches);
    CPPUNIT_TEST(testSkipInserts);
    CPPUNIT_TEST(testSkipDeletes);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
    // Keep a pointer to a single index. That index is const, so we can't mess
    // it up between test cases.
    FMDIndex const* index;
    
    // Keep a NaturalMappingScheme around.
    NaturalMappingScheme* scheme;
    
    // And a ranges bitvector of all 1s
    GenericBitVector* ranges;
    
public:
    NaturalMappingSchemeTests();
    ~NaturalMappingSchemeTests();
    
    void setUp();
    void tearDown();

    void testMap();
    void testMapWithMask();
    void testSkipMismatches();
    void testSkipInserts();
    void testSkipDeletes();
    
};

#endif
