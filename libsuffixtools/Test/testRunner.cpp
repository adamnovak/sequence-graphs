// test.cpp: Main test runner for testing the libsuffixtools library. Uses
// CppUnit. Mostly borrowed from the CPPUnit documentation at <http://cppunit.so
// urceforge.net/doc/lastest/money_example.html#sec_running_test>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <iostream>

/**
 * Main function: run all registered tests.
 */
int main(int argc, char* argv[]) {
    std::cout << "Starting" << std::endl;

    // Get the top level suite from the registry
    CppUnit::Test *suite = 
        CppUnit::TestFactoryRegistry::getRegistry().makeTest();

    // Adds the test to the list of test to run
    CppUnit::TextUi::TestRunner runner;
    runner.addTest( suite );

    // Run the tests.
    bool wasSucessful = runner.run();

    // Return error code 1 if the one of test failed.
    return wasSucessful ? 0 : 1;
}

#include <cppunit/extensions/HelperMacros.h>

/**
 * Test fixture.
 */
class BWTTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE( BWTTest );
    CPPUNIT_TEST( testConstructor );
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp();
    void tearDown();

    void testConstructor();
};

// Test fixture definition.

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( BWTTest );


void 
BWTTest::setUp()
{
}


void 
BWTTest::tearDown()
{
}


void 
BWTTest::testConstructor()
{
  CPPUNIT_FAIL( "not implemented" );
}

