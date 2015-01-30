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
 * Test looking up intervals by their left ends, from the right.
 */
void IntervalIndexTests::testLookupStartingBefore() {
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
 * Test looking up intervals by their left ends, from the left.
 */
void IntervalIndexTests::testLookupStartingAfter() {
    // First we need to define the intervals: (start, length) and a value
    std::vector<std::pair<std::pair<size_t, size_t>, std::string>> data = {
        {{2, 3}, "Ash"},
        {{3, 7}, "Bret"},
        {{10, 1}, "Corey"} 
    };
    
    // Then we make in IntervalIndex
    IntervalIndex<std::string> index(data);
    
    // Make sure the last start point exists properly
    CPPUNIT_ASSERT(index.hasStartingAfter(0));
    CPPUNIT_ASSERT(index.hasStartingAfter(1));
    CPPUNIT_ASSERT(index.hasStartingAfter(2));
    CPPUNIT_ASSERT(index.hasStartingAfter(10));
    CPPUNIT_ASSERT(!index.hasStartingAfter(11));
    
    // And we can go off the end
    CPPUNIT_ASSERT(!index.hasStartingAfter(100));
    
    // Make sure we transition between intervals where we ought to
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getStartingAfter(0).second);
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getStartingAfter(1).second);
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getStartingAfter(2).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getStartingAfter(3).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getStartingAfter(4).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getStartingAfter(5).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getStartingAfter(9).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getStartingAfter(10).second);
}

/**
 * Test looking up intervals by their right ends, from the left.
 */
void IntervalIndexTests::testLookupEndingAfter() {
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
    CPPUNIT_ASSERT(!index.hasEndingAfter(11));
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
    
    // Can't go off the end because we can't look up results we don't have.
}

/**
 * Test looking up intervals by their right ends, from the right.
 */
void IntervalIndexTests::testLookupEndingBefore() {
    // First we need to define the intervals: (start, length) and a value
    std::vector<std::pair<std::pair<size_t, size_t>, std::string>> data = {
        {{2, 3}, "Ash"},
        {{3, 7}, "Bret"},
        {{10, 1}, "Corey"} 
    };
    
    // Then we make in IntervalIndex
    IntervalIndex<std::string> index(data);
    
    // Make sure the last end point exists properly
    CPPUNIT_ASSERT(!index.hasEndingBefore(0));
    CPPUNIT_ASSERT(!index.hasEndingBefore(3));
    CPPUNIT_ASSERT(index.hasEndingBefore(4));
    CPPUNIT_ASSERT(index.hasEndingBefore(5));
    
    // And we can go off the end
    CPPUNIT_ASSERT(index.hasEndingBefore(100));
    
    // Make sure we transition between intervals where we ought to
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getEndingBefore(4).second);
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getEndingBefore(5).second);
    CPPUNIT_ASSERT_EQUAL("Ash"_s, index.getEndingBefore(8).second);
    CPPUNIT_ASSERT_EQUAL("Bret"_s, index.getEndingBefore(9).second);
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getEndingBefore(10).second);
    
    // And we can go off the end
    CPPUNIT_ASSERT_EQUAL("Corey"_s, index.getEndingBefore(100).second);
}
