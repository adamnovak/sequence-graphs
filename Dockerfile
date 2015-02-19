# This file can be used by Docker to create a container within which the
# sequence graph test code can run successfully. It can also double as
# installation instructions.

# Start with an Ubuntu that provides a decent C++11 compiler (that doesn't stub
# out std::regex, 14.04), because that's hard to change. We need gcc 4.9+
FROM ubuntu:14.10 

# I made this!
MAINTAINER Adam Novak <anovak@soe.ucsc.edu>

# Get Ubuntu Whatever with a JDK on it, and up to date packages.

# Use the commands from <https://github.com/dockerfile/java/blob/master/oracle-
# java7/Dockerfile> to install Oracle Java 7 and programmatically accept the
# license.
RUN \
    echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
    add-apt-repository -y ppa:webupd8team/java && \
    apt-get update && \
    apt-get install -y oracle-java7-installer && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/oracle-jdk7-installer

# And set JAVA_HOME    
ENV JAVA_HOME /usr/lib/jvm/java-7-oracle

# Set up to install all my packages
RUN apt-get update
RUN apt-get upgrade -y

# Get Git
RUN apt-get install -y git

# Put in GCC. We need 4.9+ for std::regex support
RUN apt-get install -y gcc g++ build-essential

# Get Boost
RUN apt-get install -y libboost-all-dev

# Get sparse hash library
RUN apt-get install -y libsparsehash-dev

# Get cppunit
RUN apt-get install -y libcppunit-dev

# Get cmake
RUN apt-get install -y cmake

# Get maven
RUN apt-get install -y maven

# Get Python 2.7
RUN apt-get install -y python2.7

# Get PCRE, which we need to build SWIG
RUN apt-get install -y libpcre3-dev

# Get SWIG 3, which we need to parse C++11. This is going to be in Ubuntu Utopic
# and higher, but we have to work with the Ubuntu that the Java images use, so
# we have to build from source.
RUN wget http://downloads.sourceforge.net/swig/swig-3.0.5.tar.gz && \
    tar -xvzf swig-3.0.5.tar.gz
RUN cd swig-3.0.5 && ./configure --prefix=/usr && make && make install

# Get SBT manually, since their HTTPS source segfaults my apt-get
RUN wget https://dl.bintray.com/sbt/debian/sbt-0.13.7.deb
RUN dpkg -i sbt-0.13.7.deb

# Get the SDSL dependency and install it under /usr so we don't need to mess
# with compiler paths.
RUN git clone https://github.com/simongog/sdsl-lite.git && \
    cd sdsl-lite && \
    ./install.sh /usr

# Clone my code and all the submodules
RUN git clone --recursive https://github.com/adamnovak/sequence-graphs.git && \
    cd sequence-graphs && \
    make

# Set everything up so we can use the container like a command
ENTRYPOINT ["sequence-graphs/createIndex/createIndex"]
CMD ["--help"]
