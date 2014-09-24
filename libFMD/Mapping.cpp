#include "Mapping.hpp"

Mapping::Mapping(): location(0, 0), is_mapped(false) {
}

Mapping::Mapping(TextPosition location, bool is_mapped, size_t context): 
    location(location), is_mapped(is_mapped), context(context) {
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
        if(context > 0) {
            o << "(+" << context << ")";
        }
    } else {
        o << "-----------------";
    }
    return o;
}
