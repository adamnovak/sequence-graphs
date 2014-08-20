# Makefile for libsuffixtools, a library of FM- and FMD-index tools for DNA
# sequences.

CXX = g++

DEBUG_FLAGS = -g -pg

HASHMAP_FLAGS = -DHAVE_EXT_HASH_MAP -DHAVE_GOOGLE_SPARSE_HASH_MAP

# We need to define CFLAGS since some of our thirdparty code is in C and needs
# -fPIC
CFLAGS = -O3 -Wall -Wextra -fPIC -pthread $(DEBUG_FLAGS) $(PARALLEL_FLAGS) \
    $(HASHMAP_FLAGS)
    
CXXFLAGS = $(CFLAGS) -std=c++11
    
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
    Alphabet.o \
    BitChar.o \
    BitVector.o \
    ClusterReader.o \
    Contig.o \
    CorrectionThresholds.o \
    DNAString.o \
    Interval.o \
    KmerDistribution.o \
    Match.o \
    MultiAlignment.o \
    MultiOverlap.o \
    Pileup.o \
    PrimerScreen.o \
    Quality.o \
    QualityTable.o \
    QualityVector.o \
    ReadInfoTable.o \
    ReadTable.o \
    SGAStats.o \
    SeqCoord.o \
    SeqReader.o \
    StdAlnTools.o \
    Util.o \
    VCFUtil.o \
    bucketSort.o \
    Thirdparty/bcr.o \
    Thirdparty/stdaln.o \
    Thirdparty/gzstream.o
    
# What do we need for our test runner binary?
TEST_OBJS = Test/testRunner.o \
	Test/suffixArrayTests.o \
	Test/sampledSuffixArrayTests.o \
	Test/mkqsTests.o

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
	
# We can automagically get header dependencies.
dependencies.mk: *.cpp *.h Thirdparty/*.C Thirdparty/*.c Thirdparty/*.h Test/*.cpp Test/*.h
	g++ $(CXXFLAGS) -MM *.cpp Thirdparty/*.C Thirdparty/*.c Test/*.cpp > dependencies.mk
	
# Include auto-generated dependencies.
include dependencies.mk
