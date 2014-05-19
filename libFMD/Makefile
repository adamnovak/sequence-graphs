CXX=g++

# What are our generic objects?
OBJS=FMDIndex.o FMDIndexBuilder.o util.o FMDIndexIterator.o Mapping.o \
	FMDPosition.o BitBuffer.o BitVector.o NibbleVector.o
    
# What do we need for our test runner binary?
TEST_OBJS=Test/TestRunner.o Test/BWTTests.o Test/FMDIndexBuilderTests.o \
    Test/FMDIndexTests.o

# What projects do we depend on? We have rules for each of these.
DEPS=libsuffixtools

# Specify all the libs to link with.
LDLIBS += ../libsuffixtools/libsuffixtools.a -lboost_filesystem

LDFLAGS += -L../libsuffixtools

# We need to link against cppunit when we run our tests.
TEST_LIBS = -lcppunit

CXXFLAGS += -O3 -fPIC -g -I../deps -I../libsuffixtools

# Stop deleting intermediate files I might need to use in the final program!
.SECONDARY:

default: all

all: libfmd.a libfmd.so

libfmd.a: $(OBJS)
	ar rcs $@ $^
	
libfmd.so: $(OBJS) $(SWIG_OBJS)
	$(CXX) $(LDFLAGS) -shared -o $@ $^
	
# libsuffixtools dependency
libsuffixtools:
	$(MAKE) -C ../libsuffixtools

clean:
	rm -Rf $(TEST_OBJS) $(OBJS) testRunner libfmd.a libfmd.so
	
test: check

check: testRunner
	./testRunner
	
testRunner: $(TEST_OBJS) $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(TEST_OBJS) $(OBJS) $(LDLIBS) \
	$(TEST_LIBS)
	
# We can automagically get header dependencies.
dependencies.mk:
	g++ $(CXXFLAGS) -MM *.cpp Test/*.cpp > dependencies.mk
	
# Include auto-generated dependencies.
include dependencies.mk
    


