# Makefile for libsuffixtools, a library of FM- and FMD-index tools for DNA
# sequences.

CXX = g++

DEBUG_FLAGS = -g

PARALLEL_FLAGS = -DHAVE_OPENMP -fopenmp

HASHMAP_FLAGS = -DHAVE_EXT_HASH_MAP -DHAVE_GOOGLE_SPARSE_HASH_MAP

CXXFLAGS = -Wall -Wextra -O3 -fPIC $(DEBUG_FLAGS) $(PARALLEL_FLAGS) \
    $(HASHMAP_FLAGS) -IUtil -IThirdparty -IConcurrency

# What object files need to be built and put in the library?
OBJS = BWTAlgorithms.o \
	BWTCABauerCoxRosone.o \
	BWTCARopebwt.o \
	BWTDiskConstruction.o \
	BWTIntervalCache.o \
	BWTReaderAscii.o \
	BWTReaderBinary.o \
	BWTReader.o \
	BWTWriterAscii.o \
	BWTWriterBinary.o \
	BWTWriter.o \
	GapArray.o \
	InverseSuffixArray.o \
	Occurrence.o \
	PopulationIndex.o \
	QuickBWT.o \
	RankProcess.o \
	RLBWT.o \
	SACAInducedCopying.o \
	SampledSuffixArray.o \
	SAReader.o \
	SAWriter.o \
	SBWT.o \
	STCommon.o \
	SuffixArray.o \
	SuffixCompare.o

# What dependency libraries do we need to include?
DEPS = Thirdparty/thirdparty.a Util/util.a

# We want those all to be packed up into libraries.

libsuffixtools.a: $(OBJS)
	ar rcs $@ $(OBJS)
	
libsuffixtools.so: $(OBJS)
	$(CXX) $(LDFLAGS) -shared -o $@  $(OBJS)
	
clean:
	rm -f *.a *.so *.o
