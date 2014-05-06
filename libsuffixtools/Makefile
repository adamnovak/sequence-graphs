# Makefile for libsuffixtools, a library of FM- and FMD-index tools for DNA
# sequences.

CXX = g++

DEBUG_FLAGS = -g

PARALLEL_FLAGS = -DHAVE_OPENMP -fopenmp

HASHMAP_FLAGS = -DHAVE_EXT_HASH_MAP -DHAVE_GOOGLE_SPARSE_HASH_MAP

# We need to define CFLAGS since some of our thirdparty code is in C and needs
# -fPIC
CFLAGS = -Wall -Wextra -O3 -fPIC -pthread $(DEBUG_FLAGS) $(PARALLEL_FLAGS) \
    $(HASHMAP_FLAGS) -IUtil -IThirdparty -IConcurrency
    
CXXFLAGS = $(CFLAGS)
    
LDFLAGS = -lz
    
# We need to link against cppunit when we run our tests.
TEST_LIBS = -lcppunit
   

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
    SuffixCompare.o \
    Util/Alphabet.o \
    Util/BitChar.o \
    Util/BitVector.o \
    Util/BloomFilter.o \
    Util/ClusterReader.o \
    Util/Contig.o \
    Util/CorrectionThresholds.o \
    Util/DNAString.o \
    Util/Interval.o \
    Util/KmerDistribution.o \
    Util/Match.o \
    Util/MultiAlignment.o \
    Util/MultiOverlap.o \
    Util/Pileup.o \
    Util/PrimerScreen.o \
    Util/Quality.o \
    Util/QualityTable.o \
    Util/QualityVector.o \
    Util/ReadInfoTable.o \
    Util/ReadTable.o \
    Util/SGAStats.o \
    Util/SeqCoord.o \
    Util/SeqReader.o \
    Util/StdAlnTools.o \
    Util/Util.o \
    Util/VCFUtil.o \
    Util/bucketSort.o \
    Thirdparty/MurmurHash3.o \
    Thirdparty/bcr.o \
    Thirdparty/stdaln.o \
    Thirdparty/gzstream.o
    
# What do we need for our test runner binary?
TEST_OBJS = Test/testRunner.o Test/bwtTest.o

default: libsuffixtools.a libsuffixtools.so

test: check

check: testRunner
	./testRunner

# We want those all to be packed up into libraries.

libsuffixtools.a: $(OBJS)
	ar rcs $@ $(OBJS)
	
libsuffixtools.so: $(OBJS)
	$(CXX) $(LDFLAGS) -shared -o $@ $(OBJS)
	
testRunner: $(TEST_OBJS) libsuffixtools.a
	$(CXX) $(CXXFLAGS) -o $@ $^ libsuffixtools.a $(LDFLAGS) \
	$(TEST_LIBS)
	
clean:
	rm -f *.a *.so *.o */*.o testRunner
