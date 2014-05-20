CXX=g++

# What are our generic objects?
OBJS=FMDIndex.o FMDIndexBuilder.o util.o FMDIndexIterator.o Mapping.o \
	FMDPosition.o BitBuffer.o BitVector.o NibbleVector.o
	
# Waht are our SWIG JNI wrapper objects?
SWIG_OBJS=swigbindings_wrap.o
    
# What do we need for our test runner binary?
TEST_OBJS=Test/TestRunner.o Test/BWTTests.o Test/FMDIndexBuilderTests.o \
    Test/FMDIndexTests.o

# What projects do we depend on? We have rules for each of these.
DEPS=libsuffixtools

# Specify all the libs to link with.
LDLIBS += ../libsuffixtools/libsuffixtools.a -lboost_filesystem -lboost_system

LDFLAGS += -L../libsuffixtools

# We need to link against cppunit when we run our tests.
TEST_LIBS = -lcppunit

CXXFLAGS += -O3 -fPIC -g -I../deps -I../libsuffixtools

# What Java package should we put the SWIG bindings in? Also used as the Maven
# groupID.
JAVA_PACKAGE = org.ga4gh

# Flags to use for building SWIG code against Java. Adjust for your platform
SWIG_FLAGS = -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -fno-strict-aliasing

# Stop deleting intermediate files I might need to use in the final program!
.SECONDARY:

default: all

all: libfmd.a libfmd.so

libfmd.a: $(OBJS)
	ar rcs $@ $^
	
libfmd.so: $(OBJS) $(SWIG_OBJS)
	$(CXX) $(LDFLAGS) -shared -o $@ $^ $(LDLIBS)
	
# libsuffixtools dependency
libsuffixtools:
	$(MAKE) -C ../libsuffixtools
	
jar: libfmd.jar

libfmd.jar: libfmd.so FMDNativeLoader.java
	mkdir -p jar
	javac java/*.java FMDNativeLoader.java -d jar
	# Make the directory for the Java package
	mkdir -p jar/`echo "$(JAVA_PACKAGE)" | sed s/\\\\./\\\\//g`
	cp libfmd.so jar/`echo "$(JAVA_PACKAGE)" | sed s/\\\\./\\\\//g`/
	jar cf $@ -C jar .

# Install the jar in the Maven local repository. See
# http://maven.apache.org/guides/mini/guide-3rd-party-jars-local.html
jar-install: libfmd.jar
	mvn install:install-file -Dfile=libfmd.jar -DgroupId=$(JAVA_PACKAGE) \
	-DartifactId=libfmd -Dversion=1.0.0-SNAPSHOT -Dpackaging=jar

# SWIG C++ file generation and compilation.	
%_wrap.o: %_wrap.cxx
	$(CXX) $(INCLUDES) $(SWIG_FLAGS) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $^
	
%_wrap.cxx: %.i
	mkdir -p java
	swig -c++ -java -outdir java -package $(JAVA_PACKAGE) $(SIZE_FLAGS) $(VECTOR_FLAGS) $^

clean:
	rm -Rf *.o testRunner libfmd.a libfmd.so libfmd.jar java/ jar/ *_wrap.cxx
	
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
    


