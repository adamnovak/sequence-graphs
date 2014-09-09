#include "Fasta.hpp"

#include <cctype>
#include <stdexcept>

Fasta::Fasta(std::string filename): stream(filename.c_str()) {
    // Nothing to do; file is open.
}

bool Fasta::hasNext() {
    if(!stream.good()) {
        // Something bad happened (like an EOF)
        return false;
    }
    
    // Look at the next character
    int peeked = stream.peek();
    while(peeked != std::char_traits<char>::eof() && isspace(peeked)) {
        // Scan ahead until it's something important.
        stream.get();
        peeked = stream.peek();
    }
    
    if(!stream.good()) {
        // Something bad happened (like an EOF)
        return false;
    }
    
    if(peeked != '>') {
        // Somebody wrote something silly, like a FASTA comment
        return false;
    }
    
    
    return true;
}

std::string Fasta::getNext() {
    return getNextRecord().second;
}

std::pair<std::string, std::string> Fasta::getNextRecord() {
    // First we have to read the header
    std::string header;
    
    // Grab the '>'
    int got = stream.get();
    if(got != '>') {
        // Complain about an error in the FASTA
        throw std::runtime_error(std::string(
            "FASTA parser expected: '>', got '") + (char)got + "'");
    }
    
    if(!stream.good()) {
        // Complain if something happens to the stream after that character.
        throw std::runtime_error(
            "Ran out of data after parsing record start character");
    }
    
    // Keep track of the character we just got
    got = stream.get();
    while(stream.good() && got != '\n') {
        // Get characters until we hit the end of the first (header) line.
        header.push_back(got);
        got = stream.get();
    }
    
    // We're going to fill up this string.
    std::string sequence = "";
    
    // Keep track of the next character.
    int peeked = stream.peek();
    while(stream.good() && peeked != '>' && peeked != 
        std::char_traits<char>::eof()) {
        
        // For each character until the next header or EOF
        
        // Grab it
        got = stream.get();
        
        if(!isspace(got)) {
            // Put it in the record string, skipping whitespace
            sequence.push_back((char) got);
        }
        
        // Look at the next character.
        peeked = stream.peek();
    }
    
    // Return our built record string.
    return std::make_pair(header, sequence);
}
