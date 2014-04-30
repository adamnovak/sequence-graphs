##Dependencies

###To Install Spark 1.0.0

Make sure to set `SPARK_HADOOP_VERSION=<your Hadoop version>` and `SPARK_YARN=<true or false depending on whether you want to use the new map-reduce>`.

```
export SPARK_HADOOP_VERSION=<your Hadoop version>
export SPARK_YARN=<whether you want to use YARN or not>

git clone git@github.com:apache/incubator-spark.git
cd incubator-spark
sbt clean
sbt assembly
sbt publish-local
```

### To Install ADAM

```
git clone git@github.com:bigdatagenomics/adam.git
cd adam
```

Make sure to edit `pom.xml`, following the instructions to fill in your Hadoop version. If you are running against a Mesos cluster, it is vitally important to be correct about whether Hadoop 2.2.0 or greater is being used or not, because this determines the version of protobuf that gets used. An incorrect version of protobuf can lead to jobs segfaulting when trying to load Mesos

```
MAVEN_OPTS="-Xmx512m -XX:MaxPermSize=128m" mvn clean package install
```

### To Install Boost 1.38 or higher, including the filesystem, system and program_options libraries

This really should come from your distribution; C++ without Boost isn't a real language anymore. If you want to install it yourself, refer to the Boost documentation. If you install it in a non-standard location, be sure to set `LIBRARY_PATH`, `CPLUS_INCLUDE_PATH`, `C_INCLUDE_PATH`, `LD_LIBRARY_PATH`, and `LD_RUN_PATH` as appropriate so that your gcc knows where to find it.

### To Install Avro 1.7.6 for C++

```
wget http://mirror.cogentco.com/pub/apache/avro/avro-1.7.6/avro-src-1.7.6.tar.gz
tar -xvzf avro-src-1.7.6.tar.gz
cd avro-src-1.7.6
cd lang/c++
./build.sh test
sudo ./build.sh install
```

If you are not root, you will have to install Avro to a nonstandard location. If your standard nonstandard location is `~/.local`, you can do something like:

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

##Installation

###Cloning

```
git clone https://github.com/adamnovak/sequence-graphs.git
cd sequence-graphs
git submodule init
git submodule update
```

###Building

####Building RLCSA

Mapping into sequence graphs relies on the Run-Length-Compressed Suffix Arrays (RLCSA) library. More specifically, a fork of the library with FMD-Index support and the appropriate SWIG bindings is required. This library ships with the program, but needs to be installed into your local Maven repository, and have some tools available on your `$PATH`. To install RLCSA's SWIG bindings in your local Maven repository:

```
cd deps/rlcsa
make
make rlcsa_grep
make jar-install
```

If this doesn't work, make sure you have SWIG, Maven, and JDK 7+ installed, that your $JAVA_HOME environment variable is set, and that you are building on Linux. Building and loading natives for other platforms (OS X) will require enhancements to the current RLCSA build system.

Once this is done, make sure to add the RLCSA tools to your `$PATH` in your `.bashrc`. If you cloned into your home directory, that would be:

```
export PATH=$PATH:$HOME/sequence-graphs/deps/rlcsa
```

####Building createIndex

The tool to create sequence graph indexes for mapping to, `createIndex`, needs to be built manually since I don't know how to hook up makefiles to SBT yet. From the repository root, if you've done everything above, all you have to do is:

```
cd createIndex
make
```

The tool can then be run as `createIndex` in that directory, or as `createIndex.sh` in the repository root.

####Building the Scala Tools

All the Scala tools and libraries can be built with just one command:

```
sbt stage
```

###Testing

```
sbt test
```

###Packaging

A single JAR which includes all dependencies can be created by running:

```
sbt assembly
```

###Running command-line tools

```
./createIndex.sh <index directory name> [<fasta> [<fasta> [<fasta> ... ]]]

./mapToIndex.sh [--cluster <Spark cluster URL>] [--repeat <times to repeat mapping>] <index directory name> <literal DNA sequence>
```
