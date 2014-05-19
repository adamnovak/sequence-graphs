#ifndef CREATEINDEX_UTIL_HPP
#define CREATEINDEX_UTIL_HPP

#include <string>

/**
 * File of useful utility functions.
 */
 
/**
 * How many bases are there? Only 4 in SGA world; N isn't allowed.
 */
static const size_t NUM_BASES = 4;

/**
 * This holds the bases in alphabetcal order by reverse complement. This order
 * is needed when doing the iterative scoping out of the reverse complement
 * intervals in the extension procedure, and there we need to go through them in
 * this order. See <http://stackoverflow.com/q/2312860/402891>
 */
extern const std::string BASES;

/**
 * This holds the bases in alphabetical order, which is the same as their sort
 * order in the BWT, and is the order that FMDIterator needs to go through them
 * in.
 */
extern const std::string ALPHABETICAL_BASES;

/**
 * Return true if a character is a valid DNA base, and false otherwise.
 */
inline bool isBase(char input) {
    for(std::string::const_iterator i = BASES.begin(); i != BASES.end(); ++i) {
        if(input == *i) {
            return true;
        }
    }
    return false;
}

// We provide a complement(char) and reverseComplement(std::string&) from
// libsuffixtools.
#include "Util.h"

/**
 * Make a new temporary directory securely, in the appropriate directory
 * (obeying $TMPDIR if set), and return its path.
 */
std::string make_tempdir();

#endif
