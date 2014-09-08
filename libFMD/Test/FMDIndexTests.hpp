#ifndef FMDINDEXBUILDERTESTS_HPP
#define FMDINDEXBUILDERTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../FMDIndex.hpp"

/**
 * Tests for the FMDIndex.
 */
class FMDIndexTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(FMDIndexTests);
    CPPUNIT_TEST(testDump);
    CPPUNIT_TEST(testDisplay);
    CPPUNIT_TEST(testMetadata);
    CPPUNIT_TEST(testLF);
    CPPUNIT_TEST(testSearch);
    CPPUNIT_TEST(testLocate);
    CPPUNIT_TEST(testIterate);
    CPPUNIT_TEST(testDisambiguate);
    CPPUNIT_TEST(testMap);
    CPPUNIT_TEST(testContextLimit);
    CPPUNIT_TEST(testAddContext);
    CPPUNIT_TEST(testMultContext);
    CPPUNIT_TEST(testLCP);
    CPPUNIT_TEST(testRetract);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
    // Keep a pointer to a single index. That index is const, so we can't mess
    // it up between test cases.
    FMDIndex const* index;
    
public:
    FMDIndexTests();
    ~FMDIndexTests();
    
    void setUp();
    void tearDown();

    void testMetadata();
    void testLF();
    void testDump();
    void testDisplay();
    void testSearch();
    void testLocate();
    void testIterate();
    void testDisambiguate();
    void testMap();
    void testContextLimit();
    void testAddContext();
    void testMultContext();
    void testLCP();
    void testRetract();
    
};

#endif
