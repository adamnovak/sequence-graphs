// Test IntervalIndex objects.

#include "../IntervalIndex.hpp"

#include "IntervalIndexTests.hpp"

#include <vector>
#include <string>

// We want string literals! Unfortunately they're in c++14. So we make our own.
std::string operator""_s(const char* string, size_t len) {
    return std::string(string, len);
}

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
        {{2, 3}, "Ash"},
        {{3, 7}, "Bret"},
        {{10, 1}, "Corey"} 
    };
    
    // Then we make in IntervalIndex
    IntervalIndex<std::string> index(data);
}

/**
 * Test looking up intervals by their left ends.
 */
void IntervalIndexTests::testLookupBefore() {
    // First we need to define the intervals: (start, length) and a value
    std::vector<std::pair<std::pair<size_t, size_t>, std::string>> data = {
        {{2, 3}, "Ash"},
        {{3, 7}, "Bret"},
        {{10, 1}, "Corey"} 
    };
    
    // Then we make in IntervalIndex
    IntervalIndex<std::string> index(data);
    
    // Make sure the first start point exists properly
    CPPUNIT_ASSERT(!index.hasStartingBefore(0));
    CPPUNIT_ASSERT(!index.hasStartingBefore(1));
    CPPUNIT_ASSERT(index.hasStartingBefore(2));
    CPPUNIT_ASSERT(index.hasStartingBefore(3));
    CPPUNIT_ASSERT(index.hasStartingBefore(4));
    
    // And we can go off the end
    CPPUNIT_ASSERT(index.hasStartingBefore(100));
    
    // Make sure we transition between intervals where we ought to
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getStartingBefore(2).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getStartingBefore(3).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getStartingBefore(4).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getStartingBefore(5).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getStartingBefore(9).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getStartingBefore(10).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getStartingBefore(11).second);
    
    // And we can go off the end
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getStartingBefore(100).second);
}

/**
 * Test looking up intervals by their right ends.
 */
void IntervalIndexTests::testLookupAfter() {
    // First we need to define the intervals: (start, length) and a value
    std::vector<std::pair<std::pair<size_t, size_t>, std::string>> data = {
        {{2, 3}, "Ash"},
        {{3, 7}, "Bret"},
        {{10, 1}, "Corey"} 
    };
    
    // Then we make in IntervalIndex
    IntervalIndex<std::string> index(data);
    
    // Make sure the last end point exists properly
    CPPUNIT_ASSERT(index.hasEndingAfter(0));
    CPPUNIT_ASSERT(index.hasEndingAfter(10));
    CPPUNIT_ASSERT(index.hasEndingAfter(11));
    CPPUNIT_ASSERT(!index.hasEndingAfter(12));
    
    
    // And we can go off the end
    CPPUNIT_ASSERT(!index.hasEndingAfter(100));
    
    // Make sure we transition between intervals where we ought to
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getEndingAfter(0).second);
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getEndingAfter(1).second);
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getEndingAfter(2).second);
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getEndingAfter(3).second);
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getEndingAfter(4).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getEndingAfter(5).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getEndingAfter(8).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getEndingAfter(9).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getEndingAfter(10).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getEndingAfter(11).second);
    
    // Can't go off the end because we can't look up results we don't have.
}
