#ifndef MARKOVMODELTESTS_HPP
#define MARKOVMODELTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "../MarkovModel.hpp"

/**
 * Tests for the Markov model.
 */
class MarkovModelTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(MarkovModelTests);
    CPPUNIT_TEST(testCodingCost);
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the model to test with.
    static const std::string filename;
    
public:
    MarkovModelTests();
    ~MarkovModelTests();
    
    void setUp();
    void tearDown();

    void testCodingCost();
    
};

#endif
