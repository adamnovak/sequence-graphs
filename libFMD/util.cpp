#include <algorithm>
#include <cstdlib>
#include <stdexcept>

#include "util.hpp"

// Bases ordered alphabetically by reverse complement
const std::string BASES = "TGCA";
// Bases ordered alphabetically by themselves.
const std::string ALPHABETICAL_BASES = "ACGT";

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
