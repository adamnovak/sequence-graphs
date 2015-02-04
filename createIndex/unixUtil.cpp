#include "unixUtil.hpp"


#include <Log.hpp>

#include <csignal>
#include <sys/resource.h>
#include <execinfo.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

// Shouls we try to demangle C++ function names in stack traces?
#define DEMANGLE_NAMES

#ifdef DEMANGLE_NAMES
// Needed for automatic name demangling, but not all that portable
#include <cxxabi.h>
#endif

void logMemory() {
    
    // We have to interrogate /proc/self/status
    
    std::ifstream statusStream("/proc/self/status");
    
    std::string line;
    while(std::getline(statusStream, line)) {
        if(line.size() >= 2 && line[0] == 'V' && line[1] == 'm') {
            // This is a virtual memory line. Log it.
            Log::output() << line << std::endl;
        }
    }
    
    

}

void exitOnSignal(int signalNumber) {
    // Log the signal.
    Log::critical() << "Exiting on signal " << signalNumber << std::endl;
    
    // Call exit and report the signal.
    exit(signalNumber);
}

void stacktraceOnSignal(int signalNumber) {
    // How many frames can we handle?
    const size_t MAX_FRAMES = 100;
    
    // This holds the stack frames
    void *frames[MAX_FRAMES];
    
    // And this holds how many there actually are, which comes out of the
    // function that gets the frames.
    size_t framesUsed = backtrace(frames, MAX_FRAMES);
    
    // Log the signal.
    Log::critical() << "Critical signal " << signalNumber << 
        ", stacktracing." << std::endl;
        
    char** traceMessages = backtrace_symbols(frames, framesUsed);
    
    for(size_t i = 0; i < framesUsed; i++) {
        
        // Start off the stack trace message, and save the stream.
        auto out = Log::critical() << "Frame " << i << ": ";
    
        #ifdef DEMANGLE_NAMES
            // Demangle the name as seen at
            //<http://panthema.net/2008/0901-stacktrace-demangled/>
            
            // Name is module(function+offset) [address] in standard format.
            // For example: 
            // ../createIndex/createIndex(_Z12make_tempdirv+0x1a4) [0x46e8f4]
            // We need to find the start and end of the function part. Sometimes
            // the parens can just be empty, so we need to handle that too.
            
            // Unfortunately, the filename on the left can contain whatever
            // characters it wants. We need to parse form the right, so we need
            // the string length.
            size_t frameLength = strlen(traceMessages[i]);
            
            // Where is the close paren, reading from the right?
            size_t closeParen = 0;
            // Where is the plus, reading from the right? If there is no plus in
            // the parens, set to 0.
            size_t plus = 0;
            // Where is the open paren, reading from right to left?
            size_t openParen = 0;
            
            for(size_t j = frameLength - 1; j != (size_t) -1; j--) {
                // Scan from right to left.
                
                if(closeParen == 0 && traceMessages[i][j] == ')') {
                    // We found the rightmost close paren
                    closeParen = j;
                } else if(j < closeParen && plus == 0 &&
                    traceMessages[i][j] == '+') {
                    
                    // We found the + to the left of the close paren.
                    plus = j;
                } else if(j < closeParen && openParen == 0 &&
                    traceMessages[i][j] == '(') {
                        
                    // We found the open paren to the left of the close paren.
                    openParen = j;
                    
                    // We're done parsing.
                    break;
                }
                
                // TODO: also pull out the address?
            }
            
            if(openParen == 0 || closeParen == 0 || plus == 0) {
                // We couldn't pull out a name and address. Either we have a
                // nonstandard format or we have empty parens.
                
                // Just use the default trace message
                out << traceMessages[i];
            } else {
                // We did parse out stuff!
                
                for(size_t j = 0; j <= openParen; j++) {
                    // Write everything up through the open paren
                    out << traceMessages[i][j];
                }
                
                // Put a null at the end of the function name
                traceMessages[i][plus] = '\0';
                
                // Make a place for the demangling function to save its status
                int status;
                
                // Do the demangling
                char* demangled = abi::__cxa_demangle(
                    &traceMessages[i][openParen + 1], NULL, NULL, &status);
                
                if(status != 0) {
                    // If we couldn't demangle the name, just use the mangled
                    // name. It's still null-terminated.
                    demangled = &traceMessages[i][openParen + 1];
                }
                
                // Write the (probably) demangled name, a "+", and the rest of
                // the message.
                out << demangled << '+' << &traceMessages[i][plus + 1];
                
                if(status == 0) {
                    // We got a demangled name we need to clean up.
                    free(demangled);
                }
            }
            
        #else
            
            // We're not doing demangling. Just write the message.
            out << traceMessages[i];
        
        #endif
    
        // End the line after the message.
        out << std::endl;
    }
    
    // Free our stacktrace memory.
    free(traceMessages);
    
    // Dump a memory usage log too.
    logMemory();
    
    exit(signalNumber);
}

