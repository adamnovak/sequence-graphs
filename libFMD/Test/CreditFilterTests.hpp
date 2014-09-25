#ifndef CREDITFILTERTESTS_HPP
#define CREDITFILTERTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../CreditFilter.hpp"

/**
 * Tests for the CreditFilter.
 */
class CreditFilterTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(CreditFilterTests);
    CPPUNIT_TEST(testApply);
    CPPUNIT_TEST(testDisagreement);
    CPPUNIT_TEST(testDistance);
    CPPUNIT_TEST(testConflictingCredit);
    CPPUNIT_TEST(testConflictingCreditOneSideOnly);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
    // Keep a pointer to a single index. That index is const, so we can't mess
    // it up between test cases.
    FMDIndex const* index;
    
public:
    CreditFilterTests();
    ~CreditFilterTests();
    
    void setUp();
    void tearDown();

    void testApply();
    void testDisagreement();
    void testDistance();
    void testConflictingCredit();
    void testConflictingCreditOneSideOnly();
    
};

#endif
