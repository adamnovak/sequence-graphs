#include "Fasta.hpp"

#include <Log.hpp>

#include <cctype>

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
    Log::info() << "Peeked: " << peeked << "(" << (char)peeked << ")" <<
        std::endl;
    while(peeked != std::char_traits<char>::eof() && isspace(peeked)) {
        // Scan ahead until it's something important.
        stream.get();
        peeked = stream.peek();
        Log::info() << "Peeked: " << peeked << "(" << (char)peeked << ")" <<
            std::endl;
    }
    
    if(!stream.good()) {
        // Something bad happened (like an EOF)
        return false;
    }
    
    if(peeked != '>') {
        // Somebody wrote something silly, like a FASTA comment
        return false;
    }
    
    Log::info() << "Next exists!" << std::endl;
    
    return true;
}

std::string Fasta::getNext() {
    // Keep track of the character we just got
    int got = stream.get();
    Log::info() << "Got: " << got << "(" << (char)got << ")" <<
        std::endl;
    while(stream.good() && got != '\n') {
        // Get characters until we hit the end of the first (header) line.
        got = stream.get();
        Log::info() << "Got: " << got << "(" << (char)got << ")" <<
            std::endl; 
    }
    
    // We're going to fill up this string.
    std::string toReturn = "";
    
    // Keep track of the next character.
    int peeked = stream.peek();
    Log::info() << "Peeked: " << peeked << "(" << (char)peeked << ")" <<
        std::endl;
    while(stream.good() && peeked != '>' && peeked != 
        std::char_traits<char>::eof()) {
        
        // For each character until the next header or EOF
        
        // Grab it
        got = stream.get();
        
        if(!isspace(got)) {
            // Put it in the record string, skipping whitespace
            toReturn.push_back((char) got);
        }
        
        // Look at the next character.
        peeked = stream.peek();
        Log::info() << "Peeked: " << peeked << "(" << (char)peeked << ")" <<
            std::endl;
    }
    
    Log::info() << "Found " << toReturn << std::endl;
    
    // Return our built record string.
    return toReturn;
    
}
