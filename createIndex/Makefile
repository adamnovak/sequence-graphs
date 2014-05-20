CXX=g++

# What are our generic objects?
OBJS=
    
# What objects do we need for our createIndex binary?
CREATEINDEX_OBJS=createIndex.o

# What projects do we depend on? We have rules for each of these.
DEPS=pinchesAndCacti sonLib libsuffixtools libfmd

# Specify all the libs to link with.
LDLIBS += ../deps/sonLib/lib/stPinchesAndCacti.a ../deps/sonLib/lib/sonLib.a \
../libFMD/libfmd.a ../libsuffixtools/libsuffixtools.a \
-lboost_filesystem -lavrocpp

# We need to make sure all our dependency header files are where our other
# dependency includes want them to be. pinchesAndCacti just includes "sonLib.h",
# so we need to explicitly point at its include directory. And similarly we need
# to grab all the internal include directories from libsuffixtools.
CXXFLAGS += -O3 -g -I../deps -I../deps/pinchesAndCacti/inc \
-I../deps/sonLib/C/inc -I../libFMD -I../libsuffixtools

# What do we use to explode Avro IDLs into Avro schemas?
IDL2SCHEMATA=java -jar avro-tools-1.7.6.jar idl2schemata

# Stop deleting intermediate files I might need to use in the final program!
.SECONDARY:

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
libfmd:
	$(MAKE) -C ../libFMD

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
	

createIndex: $(CREATEINDEX_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(CREATEINDEX_OBJS) $(OBJS) $(LDLIBS)
	
clean:
	rm -Rf $(CREATEINDEX_OBJS) $(OBJS) createIndex \
	schemas
	
# We can automagically get header dependencies.
dependencies.mk:
	g++ $(CXXFLAGS) -MM *.cpp > dependencies.mk
	
# Include auto-generated dependencies.
include dependencies.mk
    


