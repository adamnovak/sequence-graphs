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

# What do we use to explode Avro IDLs into Avro schemas? You need to edit this
# for your system if you want to change the Avro schemas.
IDL2SCHEMATA=java -jar \
$(HOME)/build/avro-src-1.7.6/lang/java/tools/target/avro-tools-1.7.6.jar \
idl2schemata

all: createIndex

# Build the RLCSA dependency.
rlcsa:
	$(MAKE) -C ../deps/rlcsa
	
# pinchesAndCacti dependency
pinchesAndCacti:
	$(MAKE) -C ../deps/pinchesAndCacti

# Build Avro schema JSON files. IDL2SCHEMATA must be correct.
schemas/%.avsc: ../src/main/avro/reference.avdl
	mkdir -p schemas
	$(IDL2SCHEMATA) $< schemas
	
# Build Avro code-generated headers; avrogencpp must be on the PATH. And in
# order for the headers to build, the Avro headers must be where the compiler
# expects them or on CPLUS_INCLUDE_PATH.
schemas/%.hpp: schemas/%.avsc
	avrogencpp -i $< -o $@

FMDIndexBuilder.o: FMDIndexBuilder.cpp FMDIndexBuilder.hpp util.hpp kseq.h

util.o: util.cpp util.hpp

createIndex.o: createIndex.cpp schemas/Side.hpp

createIndex: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf $(OBJS) createIndex schemas
    


