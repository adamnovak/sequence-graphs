// test.cpp: Main test runner for testing the libsuffixtools library. Uses
// CppUnit. Mostly borrowed from the CPPUnit documentation at <http://cppunit.so
// urceforge.net/doc/lastest/money_example.html#sec_running_test>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <iostream>
#include <string>

/**
 * Main function: run all registered tests.
 */
int main(int argc, char** argv) {
    // Dump our arguments to avoid unused warnings, and so we can see them.
    std::cout << "Command line: ";
    for(int i = 0; i < argc; i++) {
        std::cout << argv[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Running unit tests..." << std::endl;

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

