#ifndef CREATEINDEX_UTIL_HPP
#define CREATEINDEX_UTIL_HPP

#include <string>

/**
 * File of useful utility functions.
 */
 
/**
 * Function to take the reverse complement of a DNA sequence, represented as a
 * string. The string must conatin only upper-case A, C, G, T, and N.
 */
std::string reverse_complement(const std::string& sequence);

#endif
