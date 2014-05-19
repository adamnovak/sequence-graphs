CXX=g++

# What are our generic objects?
OBJS=FMDIndex.o FMDIndexBuilder.o util.o FMDIndexIterator.o Mapping.o \
	FMDPosition.o
    
# What objects do we need for our createIndex binary?
CREATEINDEX_OBJS=createIndex.o
    
# What do we need for our test runner binary?
TEST_OBJS=test/TestRunner.o test/BWTTests.o test/FMDIndexBuilderTests.o \
    test/FMDIndexTests.o

# What projects do we depend on? We have rules for each of these.
DEPS=rlcsa pinchesAndCacti sonLib libsuffixtools

# Specify all the libs to link with.
LDLIBS += ../deps/sonLib/lib/stPinchesAndCacti.a ../deps/sonLib/lib/sonLib.a \
../libsuffixtools/libsuffixtools.a -lboost_filesystem -lrlcsa -lavrocpp

# We need -fopenmp to link against rlcsa.
LDFLAGS += -L../deps/rlcsa -fopenmp -L../deps/sonLib/lib -L../libsuffixtools

# We need to link against cppunit when we run our tests.
TEST_LIBS = -lcppunit

# We need to make sure all our dependency header files are where our other
# dependency includes want them to be. pinchesAndCacti just includes "sonLib.h",
# so we need to explicitly point at its include directory. And similarly we need
# to grab all the internal include directories from libsuffixtools.
CXXFLAGS += -O3 -g -I../deps -I../deps/pinchesAndCacti/inc \
-I../deps/sonLib/C/inc -DMASSIVE_DATA_RLCSA -I../libsuffixtools \
-I../libsuffixtools/Concurrency -I../libsuffixtools/Thirdparty \
-I../libsuffixtools/Util

# What do we use to explode Avro IDLs into Avro schemas?
IDL2SCHEMATA=java -jar avro-tools-1.7.6.jar idl2schemata

# Stop deleting intermediate files I might need to use in the final program!
.SECONDARY:

all: createIndex

# Build the RLCSA dependency.
rlcsa:
	$(MAKE) -C ../deps/rlcsa
	
# pinchesAndCacti dependency
pinchesAndCacti: sonLib
	$(MAKE) -C ../deps/pinchesAndCacti

# transitive sonLib dependency
sonLib:
	$(MAKE) -C ../deps/sonLib
	
# libsuffixtools dependency
libsuffixtools:
	$(MAKE) -C ../libsuffixtools

# Build Avro schema JSON files. IDL2SCHEMATA must be correct.
schemas/%.avsc: ../src/main/avro/reference.avdl
	mkdir -p schemas
	$(IDL2SCHEMATA) $< schemas
	
# Build Avro code-generated headers; avrogencpp must be on the PATH. And in
# order for the headers to build, the Avro headers must be where the compiler
# expects them or on CPLUS_INCLUDE_PATH.
# 
# We need to work around the fact that avrogencpp deletes the schema, yet we
# need the schema text in order to write Avro-format files.
schemas/%.hpp: schemas/%.avsc
	cp $< $<.bak
	avrogencpp -i $< -o $@
	mv $<.bak $<

# Package a schema into a char constant so we can actually use it in the code.	
schemas/%_schema.hpp: schemas/%.avsc
	cp $< $(*F)_schema
	xxd -i $(*F)_schema $@
	rm $(*F)_schema
	

FMDIndexBuilder.o: FMDIndexBuilder.cpp FMDIndexBuilder.hpp util.hpp kseq.h

FMDIndex.o: FMDIndex.cpp FMDIndex.hpp

util.o: util.cpp util.hpp

createIndex.o: createIndex.cpp schemas/Side.hpp schemas/Side_schema.hpp

createIndex: $(CREATEINDEX_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(CREATEINDEX_OBJS) $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf $(CREATEINDEX_OBJS) $(TEST_OBJS) $(OBJS) createIndex testRunner \
	schemas
	
test: check

check: testRunner
	./testRunner
	
testRunner: $(TEST_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(TEST_OBJS) $(OBJS) $(LDLIBS) \
	$(TEST_LIBS)
    


