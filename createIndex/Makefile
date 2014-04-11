CXX=g++

OBJS= createIndex.o FMDIndexBuilder.o util.o

LDLIBS += -lboost_filesystem

all: createIndex

FMDIndexBuilder.o: FMDIndexBuilder.cpp FMDIndexBuilder.hpp util.hpp kseq.h

util.o: util.cpp util.hpp

createIndex.o: createIndex.cpp

createIndex: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)
    


