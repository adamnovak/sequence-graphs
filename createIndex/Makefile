CXX=g++

# What are our generic objects?
OBJS=
    
# What objects do we need for our createIndex binary?
CREATEINDEX_OBJS=createIndex.o MergeApplier.o MergeScheme.o \
OverlapMergeScheme.o MappingMergeScheme.o

# What projects do we depend on? We have rules for each of these.
DEPS=pinchesAndCacti sonLib libsuffixtools libfmd

# Specify all the libs to link with.
LDLIBS += ../deps/sonLib/lib/stPinchesAndCacti.a ../deps/sonLib/lib/sonLib.a \
../libFMD/libfmd.a ../libsuffixtools/libsuffixtools.a \
-lboost_filesystem -lboost_program_options -lpthread -lz

# We need to make sure all our dependency header files are where our other
# dependency includes want them to be. pinchesAndCacti just includes "sonLib.h",
# so we need to explicitly point at its include directory. And similarly we need
# to grab all the internal include directories from libsuffixtools.
CXXFLAGS += -std=c++11 -O3 -g -pg -I../deps -I../deps/pinchesAndCacti/inc \
-I../deps/sonLib/C/inc -I../libFMD -I../libsuffixtools

# Stop deleting intermediate files I might need to use in the final program!
.SECONDARY:

# Re-do things every time
.PHONY: clean $(DEPS)

all: createIndex

# pinchesAndCacti dependency
pinchesAndCacti: sonLib
	$(MAKE) -C ../deps/pinchesAndCacti

# transitive sonLib dependency
sonLib:
	$(MAKE) -C ../deps/sonLib
	
# libsuffixtools dependency
libsuffixtools:
	$(MAKE) -C ../libsuffixtools
	
# libFMD dependency
libfmd: libsuffixtools
	$(MAKE) -C ../libFMD

createIndex: $(CREATEINDEX_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(CREATEINDEX_OBJS) $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf *.o createIndex
	
# We can automagically get header dependencies.
dependencies.mk: *.cpp *.hpp
	g++ $(CXXFLAGS) -MM *.cpp > dependencies.mk
	
# Include auto-generated dependencies.
include dependencies.mk
    


