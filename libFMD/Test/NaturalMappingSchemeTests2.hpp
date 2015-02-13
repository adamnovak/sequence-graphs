#ifndef NATURALMAPPINGSCHEMETESTS2_HPP
#define NATURALMAPPINGSCHEMETESTS2_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../NaturalMappingScheme.hpp"

/**
 * Test for the natural mapping scheme, using sequences we can merge.
 */
class NaturalMappingSchemeTests2 : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(NaturalMappingSchemeTests2);
    CPPUNIT_TEST(testMap);
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
    NaturalMappingSchemeTests2();
    ~NaturalMappingSchemeTests2();
    
    void setUp();
    void tearDown();

    void testMap();
    
};

#endif
