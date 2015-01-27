#ifndef LRMAPPINGSCHEMETESTS2_HPP
#define LRMAPPINGSCHEMETESTS2_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../LRMappingScheme.hpp"

/**
 * Mismatch Test for the LR MappingScheme, against a different FASTA. TODO:
 * unify somehow so onbe test suite can work on multiple FASTAs, possibly by not
 * using the fixture setup methods, or by setting up two indexes, or one index
 * of both files.
 */
class LRMappingSchemeTests2 : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(LRMappingSchemeTests2);
    CPPUNIT_TEST(testMismatch);
    CPPUNIT_TEST(testMapOnMismatch);
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
    LRMappingSchemeTests2();
    ~LRMappingSchemeTests2();
    
    void setUp();
    void tearDown();

    void testMismatch();
    void testMapOnMismatch();
    
};

#endif
