#include "Mapping.hpp"

Mapping::Mapping(): location(0, 0), is_mapped(false), context(0) {
}

Mapping::Mapping(TextPosition location, size_t context): 
    location(location), is_mapped(true), context(context) {
}

bool Mapping::operator==(const Mapping& other) const {
    return location == other.location && is_mapped == other.is_mapped && 
        context == other.context;
}

bool Mapping::operator!=(const Mapping& other) const {
    return location != other.location || is_mapped != other.is_mapped ||
        context != other.context;
}

std::ostream& operator<< (std::ostream& o, Mapping const& mapping) {
    if(mapping.is_mapped) {
        o << "Text " << mapping.location.getText() << " offset " <<
        mapping.location.getOffset();
        if(mapping.context > 0) {
            o << "(+" << mapping.context << ")";
        }
    } else {
        o << "-----------------";
    }
    return o;
}
