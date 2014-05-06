#ifndef BWTTEST_HPP
#define BWTTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

/**
 * Tests for BWT generation.
 */
class BWTTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE( BWTTest );
    CPPUNIT_TEST( testConstruction );
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;

public:
    void setUp();
    void tearDown();

    void testConstruction();
};

#endif
