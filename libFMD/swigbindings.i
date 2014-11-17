// swigbindings.i: SWIG wrapper for the FMDIndex class.

// Name the module something Java-y that isn't going to clash with the FMD
// class.
%module FMDUtil

// Set up STL compatibility
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"

// Set up int types (int64_t)
%include "stdint.i"

// Set up exception not-killing-the-process-ness. See
// <http://www.swig.org/Doc1.3/Library.html#Library_nn17>
%include "exception.i"

// Also bring in the typemaps library
%include "typemaps.i"

// Java can't handle these operator names.
%rename(operatorLeftShift) operator<<;
%rename(operatorEquals) operator==;
%rename(operatorNotEquals) operator!=;
%rename(operatorIncrement) operator++;

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// Note that build() on an FMDIndexBuilder produces a new object, which we
// should delete if we don't need it.
%newobject FMDIndexBuilder::build;

%include "CSA/BitVectorBase.hpp"
%include "CSA/BitVector.hpp"
%{
  #include "CSA/BitVectorBase.hpp"
%}
%{
  #include "CSA/BitVector.hpp"
%}

// We need to use the inner vector iterator classes to look at vectors. Give a
// partial definition under a new name.
class BitVectorIterator
{
public:
  explicit BitVectorIterator(const BitVector& par);
  ~BitVectorIterator();

  size_t rank(size_t value, bool at_least = false);

  size_t select(size_t index);
  size_t selectNext();

  pair_type valueBefore(size_t value);
  pair_type valueAfter(size_t value);
  pair_type nextValue();

  pair_type selectRun(size_t index, size_t max_length);
  pair_type selectNextRun(size_t max_length);

  bool isSet(size_t value);

  size_t countRuns();
};

// If we don't have this, it somehow manages to think that BitVector doesn't
// mean CSA::BitVector;
using CSA::BitVector;
using CSA::BitVectorEncoder;


%include "Mapping.hpp"
%include "TextPosition.hpp"
%include "MismatchResultSet.hpp"
%{
  #include "FMDIndex.hpp"
%}
%include "FMDIndex.hpp"
%include "FMDIndexIterator.hpp"
%{
  #include "FMDIndexBuilder.hpp"
%}
%include "FMDIndexBuilder.hpp"

%{
  using namespace CSA;
%}


// Since we will need to load and save range vectors to files, we need to expose
// a minimal C FILE API.
FILE* fopen(char* filename, char* mode);
void fclose(FILE* file);

// Java needs to work with vectors of mappings coming back from the map method.
%template(MappingVector) std::vector<Mapping>; 

// Java also needs to work with vectors of int64_ts coming back from the map
// method when working on ranges.
%template(IntVector) std::vector<long long>;

%template(IntPair) std::pair<int64_t, size_t>;
%template(sizePair) std::pair<size_t, size_t>;
%template(building) std::pair<int64_t, std::pair<size_t, size_t>>;
%template(IntPairVector) std::vector<std::pair<int64_t, size_t>>;
%template(IntPairPairVector) std::vector<std::pair<int64_t, std::pair<size_t, size_t>>>;

// Whenever any of the JNI classes loads, load the native library.
%pragma(java) jniclasscode=%{
  static {
    FMDNativeLoader.load();
  }
%}

// We already worked around the inner classes thing.
#pragma SWIG nowarn=SWIGWARN_PARSE_NESTED_CLASS
