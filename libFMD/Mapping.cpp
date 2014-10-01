#include "Mapping.hpp"

Mapping::Mapping(): location(0, 0), range(-1), is_mapped(false),
    leftMinContext(0), rightMinContext(0), leftMaxContext(0), 
    rightMaxContext(0) {
    
    // Nothing to do
}

Mapping::Mapping(TextPosition location): location(location), range(-1), 
    is_mapped(true), leftMinContext(0), rightMinContext(0), leftMaxContext(0), 
    rightMaxContext(0) {
    
    // Nothing to do
}

Mapping::Mapping(TextPosition location, size_t leftContext, 
    size_t rightContext): location(location), range(-1), is_mapped(true),
    leftMinContext(leftContext), rightMinContext(rightContext),
    leftMaxContext(leftContext), rightMaxContext(rightContext) {
    
    // Nothing to do
}

Mapping::Mapping(int64_t range): location(0, 0), range(range), is_mapped(true),
    leftMinContext(0), rightMinContext(0), leftMaxContext(0), 
    rightMaxContext(0) {
    
    // Nothing to do
}

bool Mapping::operator==(const Mapping& other) const {
    return location == other.location &&
        range == other.range &&
        is_mapped == other.is_mapped && 
        leftMinContext == other.leftMinContext &&
        rightMinContext == other.rightMinContext &&
        leftMaxContext == other.leftMaxContext &&
        rightMaxContext == other.rightMaxContext;
}

bool Mapping::operator!=(const Mapping& other) const {
    return !(*this == other);
}

std::ostream& operator<< (std::ostream& o, Mapping const& mapping) {
    if(mapping.is_mapped) {
        o << "Text " << mapping.location.getText() << " offset " <<
        mapping.location.getOffset();
        if(mapping.getLeftMaxContext() > 0 || 
            mapping.getRightMaxContext() > 0) {
            
            o << "(+" << mapping.leftMaxContext << "|+" << 
            mapping.rightMaxContext << ")";
        }
    } else {
        o << "-----------------";
    }
    return o;
}
