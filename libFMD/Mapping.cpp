#include "Mapping.hpp"

Mapping::Mapping(): location(0, 0), is_mapped(false), leftContext(0), 
    rightContext(0) {
    
    // Nothing to do
}

Mapping::Mapping(TextPosition location): Mapping(location, 0, 0) {
    
    // Nothing to do
}

Mapping::Mapping(TextPosition location, size_t leftContext,
    size_t rightContext): location(location), is_mapped(true), 
    leftContext(leftContext), rightContext(rightContext) {
    
    // Nothing to do
}

bool Mapping::operator==(const Mapping& other) const {
    return location == other.location && is_mapped == other.is_mapped && 
        leftContext == other.leftContext && rightContext == other.rightContext;
}

bool Mapping::operator!=(const Mapping& other) const {
    return !(*this == other);
}

std::ostream& operator<< (std::ostream& o, Mapping const& mapping) {
    if(mapping.is_mapped) {
        o << "Text " << mapping.location.getText() << " offset " <<
        mapping.location.getOffset();
        if(mapping.getContext() > 0) {
            o << "(+" << mapping.leftContext << "|+" << mapping.rightContext <<
            ")";
        }
    } else {
        o << "-----------------";
    }
    return o;
}
