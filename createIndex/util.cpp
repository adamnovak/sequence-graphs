#include <algorithm>

#include "util.hpp"


std::string reverse_complement(const std::string& sequence) {
    // Make an empty string to fill in.
    std::string toReturn;
    
    // Make it the right size.
    toReturn.resize(sequence.size)
    
    // Do a map to write to it backwards, through a shiny new C++11 anonymous
    // function.
    std::transform(sequence.begin(), sequence.end(), toReturn.rbegin(), 
        (char c) {
            // For each character, complement it
            switch(c) {
                case 'A':
                    return 'T';
                case 'C':
                    return 'G';
                case 'G':
                    return 'C';
                case 'T':
                    return 'A';
                case 'N':
                    return 'N';
                default:
                    // Complain that's not allowed.
                    throw new std::exception("Invalid character in sequence");
            }
        });
    
    // Return the string.
    return toReturn;
}
