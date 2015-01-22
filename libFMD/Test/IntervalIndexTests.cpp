// Test IntervalIndex objects.

#include "../IntervalIndex.hpp"

#include "IntervalIndexTests.hpp"

#include <vector>
#include <string>

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( IntervalIndexTests );

void IntervalIndexTests::setUp() {
}


void IntervalIndexTests::tearDown() {
}

/**
 * Test creating an IntervalIndex.
 */
void IntervalIndexTests::testCreate() {
    
    // First we need to define the intervals: (start, length) and a value
    std::vector<std::pair<std::pair<size_t, size_t>, std::string>> data = {
        {{0, 5}, "Ash"},
        {{3, 7}, "Bret"},
        {{10, 1}, "Corey"} 
    };
    
    // Then we make in IntervalIndex
    IntervalIndex<std::string> index(data);
}
