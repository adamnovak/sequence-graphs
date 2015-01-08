CXX=g++

# What are our generic objects?
OBJS=pinchGraphUtil.o unixUtil.o indexUtil.o
    
# What objects do we need for our createIndex binary?
CREATEINDEX_OBJS=createIndex.o MergeApplier.o MergeScheme.o \
OverlapMergeScheme.o MappingMergeScheme.o adjacencyComponentUtil.o

# And for our mapReads binary?
MAPREADS_OBJS=mapReads.o

# And for our cactusMerge binary?
CACTUSMERGE_OBJS=cactusMerge.o 

# And for our evaluateMapability binary?
EVALUATEMAPABILITY_OBJS=evaluateMapability.o

# What projects do we depend on? We have rules for each of these.
DEPS=pinchesAndCacti sonLib vflib libsuffixtools libfmd

# Specify all the libs to link with.
LDLIBS += ../deps/sonLib/lib/stPinchesAndCacti.a ../deps/sonLib/lib/sonLib.a \
../libFMD/libfmd.a ../libsuffixtools/libsuffixtools.a \
../deps/vflib/lib/libvf.a \
-lboost_filesystem -lboost_program_options -lpthread -lz -lboost_system -lsdsl

# We need to make sure all our dependency header files are where our other
# dependency includes want them to be. pinchesAndCacti just includes "sonLib.h",
# so we need to explicitly point at its include directory. And similarly we need
# to grab all the internal include directories from libsuffixtools.
CXXFLAGS += -std=c++11 -O3 -g -rdynamic -I../deps \
-I../deps/pinchesAndCacti/inc -I../deps/sonLib/C/inc -I../libFMD \
-I../libsuffixtools -I../deps/vflib/include

# Stop deleting intermediate files I might need to use in the final program!
.SECONDARY:

# Re-do things every time
.PHONY: clean $(DEPS)

all: createIndex cactusMerge mapReads evaluateMapability

# pinchesAndCacti dependency
pinchesAndCacti: sonLib
	$(MAKE) -C ../deps/pinchesAndCacti

# transitive sonLib dependency
sonLib:
	$(MAKE) -C ../deps/sonLib
	
# dependency on VFLib graph matching library for analyzing rearrangement
# structures
vflib:
	$(MAKE) -C ../deps/vflib
	
# libsuffixtools dependency
libsuffixtools:
	$(MAKE) -C ../libsuffixtools
	
# libFMD dependency
libfmd: libsuffixtools
	$(MAKE) -C ../libFMD

createIndex: $(CREATEINDEX_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(CREATEINDEX_OBJS) $(OBJS) $(LDLIBS)
	
mapReads: $(MAPREADS_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(MAPREADS_OBJS) $(OBJS) $(LDLIBS)
	
cactusMerge: $(CACTUSMERGE_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(CACTUSMERGE_OBJS) $(OBJS) $(LDLIBS)
	
evaluateMapability: $(EVALUATEMAPABILITY_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(EVALUATEMAPABILITY_OBJS) $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf *.o createIndex
	
# We can automagically get header dependencies.
dependencies.mk: *.cpp *.hpp
	$(CXX) $(CXXFLAGS) -MM *.cpp > dependencies.mk
	
# Include auto-generated dependencies.
include dependencies.mk
    


