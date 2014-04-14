CXX=g++

OBJS=createIndex.o FMDIndexBuilder.o util.o

# What projects do we depend on? We have rules for each of these.
DEPS=rlcsa pinchesAndCacti

# We need -fopenmp to link against rlcsa.
LDLIBS += ../deps/sonLib/lib/stPinchesAndCacti.a ../deps/sonLib/lib/sonLib.a \
-lboost_filesystem -lrlcsa

LDFLAGS += -L../deps/rlcsa -fopenmp -L../deps/sonLib/lib

# We need to make sure all our dependency header files are where our other
# dependency includes want them to be. pinchesAndCacti just includes "sonLib.h",
# so we need to explicitly point at its include directory.
CXXFLAGS += -g -I../deps -I../deps/pinchesAndCacti/inc -I../deps/sonLib/C/inc \
-DMASSIVE_DATA_RLCSA

all: createIndex

# Build the RLCSA dependency.
rlcsa:
	$(MAKE) -C ../deps/rlcsa
	
# pinchesAndCacti dependency
pinchesAndCacti:
	$(MAKE) -C ../deps/pinchesAndCacti


FMDIndexBuilder.o: FMDIndexBuilder.cpp FMDIndexBuilder.hpp util.hpp kseq.h

util.o: util.cpp util.hpp

createIndex.o: createIndex.cpp 

createIndex: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf $(OBJS) createIndex
    


