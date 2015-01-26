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
    
    // make a runner
    CppUnit::TextUi::TestRunner runner;
    
    // Get the registry
    CppUnit::TestFactoryRegistry& registry = 
        CppUnit::TestFactoryRegistry::getRegistry();
        
    // Get the test for the registry.
    CppUnit::Test* suite = registry.makeTest();
    
    // Add the test to the list of test to run
    runner.addTest(suite);
    
    // Have all the tests succeeded so far?
    bool wasSucessful = true;
    
    if(argc > 1) {
        // Say the rest of the arguments are test names and we should run them.
        
        for(size_t i = 1; i < argc; i++) {
            // Run each test in turn.
            std::cout << "Running test " << argv[i] << "..." << std::endl;
            wasSucessful &= runner.run(std::string(argv[i]));
        }
        
    } else {
        // Run all the tests
        std::cout << "Running unit tests..." << std::endl;

        // Run the tests.
        wasSucessful &= runner.run();
    }
    
    // Return error code 1 if the one of the tests failed.
    return wasSucessful ? 0 : 1;
}
