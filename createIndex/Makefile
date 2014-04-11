CXX=g++

OBJS=createIndex.o FMDIndexBuilder.o util.o

# What projects do we depend on? We have rules for each of these.
DEPS=rlcsa

# We need -fopenmp to link against rlcsa.
LDLIBS += -lboost_filesystem -lrlcsa
LDFLAGS += -L../deps/rlcsa -fopenmp

CXXFLAGS += -g -I../deps -DMASSIVE_DATA_RLCSA

all: createIndex

# Build the RLCSA dependency.
rlcsa:
	$(MAKE) -C ../deps/rlcsa

FMDIndexBuilder.o: FMDIndexBuilder.cpp FMDIndexBuilder.hpp util.hpp kseq.h

util.o: util.cpp util.hpp

createIndex.o: createIndex.cpp 

createIndex: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf $(OBJS) createIndex
    


