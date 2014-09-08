#include "unixUtil.hpp"


#include <Log.hpp>

#include <csignal>
#include <sys/resource.h>
#include <execinfo.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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
        // Log the stack frames.
        Log::critical() << "Frame " << i << ": " << traceMessages[i] << 
            std::endl;
    }
    
    exit(signalNumber);
}

