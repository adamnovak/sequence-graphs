// Test the Markov model generation.

#include "MarkovModelTests.hpp"
#include <limits>


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( MarkovModelTests );

// Define constants
const std::string MarkovModelTests::filename = "Test/markov.model";

MarkovModelTests::MarkovModelTests() {
    
}

MarkovModelTests::~MarkovModelTests() {

}

void MarkovModelTests::setUp() {
    
}


void MarkovModelTests::tearDown() {
    
}

/**
 * Test coding cost calculation
 */
void MarkovModelTests::testCodingCost() {
    
    // Load up the model
    MarkovModel model(filename);
    
    // This is too short and thus free.
    CPPUNIT_ASSERT_EQUAL(0.0, model.encodingCost("A"));
    
    // This is impossible
    CPPUNIT_ASSERT_EQUAL(std::numeric_limits<double>::infinity(), 
        model.encodingCost("AA"));
        
    // This is guaranteed and thus free
    CPPUNIT_ASSERT_EQUAL(0.0, model.encodingCost("AB"));
    
    // This is 50% and thus 1 bit
    CPPUNIT_ASSERT_EQUAL(1.0, model.encodingCost("BA"));
    
    // This is also 50% and thus 1 bit
    CPPUNIT_ASSERT_EQUAL(1.0, model.encodingCost("BB"));
    
    // This is 50% each transition from B to anything.
    CPPUNIT_ASSERT_EQUAL(6.0, model.encodingCost("BABABABBBAB"));
    
}






