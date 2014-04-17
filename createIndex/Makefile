CXX=g++

OBJS=createIndex.o FMDIndex.o FMDIndexBuilder.o util.o

# What projects do we depend on? We have rules for each of these.
DEPS=rlcsa pinchesAndCacti


LDLIBS += ../deps/sonLib/lib/stPinchesAndCacti.a ../deps/sonLib/lib/sonLib.a \
-lboost_filesystem -lrlcsa -lavrocpp

# We need -fopenmp to link against rlcsa.
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

# Stop deleting intermediate files I might need to use in the final program!
.SECONDARY:

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

createIndex: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf $(OBJS) createIndex schemas
    


