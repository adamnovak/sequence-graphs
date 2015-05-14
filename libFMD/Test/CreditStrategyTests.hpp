#ifndef CREDITSTRATEGYTESTS_HPP
#define CREDITSTRATEGYTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../CreditStrategy.hpp"

/**
 * Test for the application of credit.
 */
class CreditStrategyTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(CreditStrategyTests);
    CPPUNIT_TEST(testCredit);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;
    
    // Also we need an index temp directory
    std::string tempDir;
    
    // Keep a pointer to a single index. That index is const, so we can't mess
    // it up between test cases.
    FMDIndex const* index;
    
    // Keep a pointer to a view of the index, so we can keep it around while the
    // CreditStrategy works on it.
    FMDIndexView* view;
    
    // And a credit applier
    CreditStrategy* credit;
    
    // And a ranges bitvector merging some positions
    GenericBitVector* ranges;
    
public:
    CreditStrategyTests();
    ~CreditStrategyTests();
    
    void setUp();
    void tearDown();

    void testCredit();
};

#endif
