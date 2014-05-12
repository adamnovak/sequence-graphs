#include "Mapping.hpp"

Mapping::Mapping(): location(0, 0), is_mapped(false)
{
}

Mapping::Mapping(pair_type location, bool is_mapped): location(location), 
  is_mapped(is_mapped)
{
}

bool
Mapping::operator==(const Mapping& other) const
{
  return location == other.location && is_mapped == other.is_mapped;
}

std::ostream&
operator<< (std::ostream& o, Mapping const& mapping)
{
  if(mapping.is_mapped)
  {
    o << "Text " << mapping.location.first << " offset " <<
      mapping.location.second;
  }
  else
  {
    o << "-----------------";
  }
  return o;
}
