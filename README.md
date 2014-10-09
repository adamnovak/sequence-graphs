#Sequence Graphs

This repository contains tools for generating and working with Reference
Structures and Reference Hierarchies.

##Dependencies

### Tool Dependencies

The code in this repository is written in C++11 and Scala, with a bit of Java.
To build it, you need:

* A C++11-supporting `g++` (GCC 4.9 works great)

* `sbt`, for building Scala projects

* The JDK, for building Java code

* Maven, for installing Scala and Java JARs into the local Maven cache.

* SWIG, for generating Java bindings for C++. You need a SWIG version that 
understands C++11 syntax (SWIG 3.0+ will work).

All these should probably come from your distribution, but everything except the
JDK can be fairly easily installed in your home directory with a `--prefix`.

###C++ Library Dependencies

The code in this repository uses a few libraries at the C++ level. These need to
be built in C++11 mode, because the C++11 ABI (in GCC) is not compatible with
the C++03 ABI in some cases.

####Install Google Sparse Hash

The indexing code borrowed from SGA uses the (Google) Sparsehash library. To
install it, you can use:

```
wget https://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz
tar -xvzf sparsehash-2.0.2.tar.gz
cd sparsehash-2.0.2
./configure --prefix=$HOME/.local
make
make install
```

Make sure that `$HOME/.local/include` is in your `CPLUS_INCLUDE_PATH`, so that
the compiler can find your newly installed headers.

####Install Boost

If the version of Boost that you can install from your distribution does not
present a C++11 ABI, you will need to install Boost yourself. (If it does, you
can just install that.)

```
wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.bz2/download
tar -xvjf boost_1_55_0.tar.bz2
cd boost_1_55_0
./bootstrap.sh --prefix=$HOME/.local
./b2 cxxflags="-std=c++11"
./b2 install
```

Make sure to add `$HOME/.local/include` to your `CPLUS_INCLUDE_PATH`, and
`$HOME/.local/lib` to your `LD_LIBRARY_PATH`.

####Install SDSL

Succinct bit vectors from the Succinct Data Structures Library (SDSL) are used. You thus need to install SDSL in order to build the code. Make sure that SDSL gets picked up by your compiler.

```
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install $HOME/.local/
```


####Install CPPUnit

The C++ code includes test cases, which you can run if you have a version of
CPPUnit installed that has been built with an appropriate compiler. You probably
don't need to build it in C++11 mode, as it doesn't use any STL in its API.

```
wget "http://downloads.sourceforge.net/project/cppunit/cppunit/1.12.1/cppunit-1.12.1.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fapps%2Fmediawiki%2Fcppunit%2Findex.php%3Ftitle%3DMain_Page&ts=1403118115&use_mirror=softlayer-dal"
tar -xvzf cppunit-1.12.1.tar.gz
cd cppunit-1.12.1
./configure --prefix=$HOME/.local
make
make install
```

Again, make sure your environment variables are set so that your compiler snd
linker will pick up your installed CPPUnit.

###Scala Library Dependencies

Almost all the Scala dependencies are automatically downloaded and installed by
`sbt`. However, some need to be built from source in order for them to be
properly configured for your environment.

####Install Apache Spark 1.0.0

Make sure to set `SPARK_HADOOP_VERSION=<your Hadoop version>` and
`SPARK_YARN=<true or false depending on whether you want to use the new map-
reduce>`.

```
export SPARK_HADOOP_VERSION=<your Hadoop version>
export SPARK_YARN=<whether you want to use YARN or not>

git clone git@github.com:apache/incubator-spark.git
cd incubator-spark
sbt clean
sbt assembly
sbt publish-local
```

##Installation

###Cloning

```
git clone https://github.com/adamnovak/sequence-graphs.git
cd sequence-graphs
git submodule init
git submodule update
```

###Building

The system consists of three C++ components (`libsuffixtools`, `libFMD`, and
`createIndex`), each of which depends on the previous one, and some Scala code
which depends on a JAR file from `libFMD`.

####Building the C++ Tools and Libraries

All the C++ code can be built from the root directory of the repository with:

```
make
```

####Satisfying the Scala->C++ Dependency

In order for the Scala code to work, you need to install the `libFMD` wrapper
JAR file into your local Maven repository, where Scala can find it:

```
cd libFMD
make jar-install
cd ..
```

####Building the Scala Tools

All the Scala tools and libraries can be built with just one command:

```
sbt stage
```

###Testing

The C++ `libsuffixtools` and `libFMD` libraries have CPPUnit tests. To invoke
them, run:

```
cd libsuffixtools
make check
cd ../libFMD
make check
```


The scala code can be tested with:

```
sbt test
```

###Packaging

A single JAR of the Scala code which includes all dependencies (including the
native `libFMD` library) can be created by running:

```
sbt assembly
```

###Running command-line tools

```
./createIndex.sh [--context <number of bases to merge on>] <index directory name> [<fasta> [<fasta> [<fasta> ... ]]]

./mapToIndex.sh [--cluster <Spark cluster URL>] [--repeat <times to repeat mapping>] <index directory name> <literal DNA sequence>
```

##Additional Procedures

Here are some other things you can install if you want to.

###Install Avro 1.7.6 for C++

This library lets you read and write Avro-format data from C++.

```
wget http://mirror.cogentco.com/pub/apache/avro/avro-1.7.6/avro-src-1.7.6.tar.gz
tar -xvzf avro-src-1.7.6.tar.gz
cd avro-src-1.7.6
cd lang/c++
export CXXFLAGS=-stc=c++11 # TODO: This may or may not work
./build.sh test
sudo ./build.sh install
```

If you are not root, you will have to install Avro to a nonstandard location. If
your standard nonstandard location is `~/.local`, you can do something like:

```
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local build
cd build && make && make install
```

And then in your .bashrc make sure you have something like:

```
export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$HOME/.local/lib:$LD_RUN_PATH
export LIBRARY_PATH=$HOME/.local/lib:$LIBRARY_PATH
export C_INCLUDE_PATH=$HOME/.local/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$HOME/.local/include:$CPLUS_INCLUDE_PATH
export PATH=$HOLE/.local/bin:$PATH
```
