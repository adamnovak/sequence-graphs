#include <algorithm>
#include <cstdlib>
#include <stdexcept>

#include "util.hpp"


/*std::string reverse_complement(const std::string& sequence) {
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
}*/

/**
 * Complement a single upper-case base, or throw an exception.
 */
inline char complement(char base) {
    switch(base) {
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
            throw std::runtime_error("Invalid character in sequence");
    }
}

std::string reverse_complement(const std::string& sequence) {
    // Make an empty string to fill in.
    std::string toReturn;
    
    // Make it the right size.
    toReturn.resize(sequence.size());
    
    // Do a map to write to it backwards, through a shiny new C++11 anonymous
    // function.
    std::transform(sequence.begin(), sequence.end(), toReturn.rbegin(),
        complement);
    
    // Return the string.
    return toReturn;
}

std::string make_tempdir() {
    // This string holds the name of a temporary directory we're using.
    std::string tempDir;
    
    if(getenv("TMPDIR") == NULL) {
        // Default temp files to /tmp
        tempDir = std::string("/tmp");
    } else {
        // Put them where they really belong
        tempDir = std::string(getenv("TMPDIR"));
    }
    
    // Turn the temp dir string into a proper mkdtemp template
    tempDir += "/tmpXXXXXX";
    // Make the directory, and modify the string in place.
    if(mkdtemp(const_cast<char*>(tempDir.c_str())) == NULL) {
        throw std::runtime_error("Could not create temporary directory " +
            tempDir);
    }
    
    // Return the modified string.
    return tempDir;
}
