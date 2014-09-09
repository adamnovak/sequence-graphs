# Instructions for building speedMap, which compares mapping speeds between
# algorithms.

CXX=g++

# What are our generic objects?
OBJS=
    
# What objects do we need for our speedMap binary?
SPEEDMAP_OBJS=speedMap.o

# What projects do we depend on? We have rules for each of these.
DEPS=libsuffixtools libfmd

# Specify all the libs to link with.
LDLIBS += ../libFMD/libfmd.a ../libsuffixtools/libsuffixtools.a \
-lboost_filesystem -lboost_program_options -lpthread -lz -lboost_system -lsdsl

# We need to make sure all our dependency header files are where our other
# dependency includes want them to be. We need to grab all the internal include
# directories from libsuffixtools.
CXXFLAGS += -std=c++11 -O3 -g -I../deps -I../libFMD -I../libsuffixtools

# Stop deleting intermediate files I might need to use in the final program!
.SECONDARY:

# Re-do things every time
.PHONY: clean $(DEPS)

all: speedMap
	
# libsuffixtools dependency
libsuffixtools:
	$(MAKE) -C ../libsuffixtools
	
# libFMD dependency
libfmd: libsuffixtools
	$(MAKE) -C ../libFMD

speedMap: $(SPEEDMAP_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SPEEDMAP_OBJS) $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf *.o speedMap
	
# We can automagically get header dependencies.
dependencies.mk: *.cpp *.hpp
	g++ $(CXXFLAGS) -MM *.cpp > dependencies.mk
	
# Include auto-generated dependencies.
include dependencies.mk
    


