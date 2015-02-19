# This file can be used by Docker to create a container within which the
# sequence graph test code can run successfully. It can also double as
# installation instructions.

# This file borrows heavvily from dockerfiles from Dockerfile.

# Start with an Ubuntu that provides a decent C++11 compiler (that doesn't stub
# out std::regex, 14.04), because that's hard to change. We need gcc 4.9+
FROM ubuntu:14.10 

# I made this!
MAINTAINER Adam Novak <anovak@soe.ucsc.edu>

# Turn on Multiverse, as in <https://github.com/dockerfile/ubuntu/blob/master/Dockerfile>
RUN sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list

# We need apt-add-repository to add the java installer repository.

# Use the commands from <https://github.com/dockerfile/java/blob/master/oracle-
# java7/Dockerfile> to install Oracle Java 7 JDK and programmatically accept the
# license. But then also install all the packages we need and clear out the
# package index.
RUN \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y software-properties-common && \
    echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
    add-apt-repository -y ppa:webupd8team/java && \
    apt-get update && \
    apt-get install -y oracle-java7-installer wget git gcc g++ build-essential \
    libboost-all-dev libsparsehash-dev libcppunit-dev  cmake maven python2.7 \
    libpcre3-dev && \
    rm -rf /var/cache/oracle-jdk7-installer && \
    rm -rf /var/lib/apt/lists/*

# And set JAVA_HOME    
ENV JAVA_HOME /usr/lib/jvm/java-7-oracle

# Package installation summary:
# Get Git
# Put in GCC. We need 4.9+ for std::regex support
# Get Boost
# Get sparse hash library
# Get cppunit
# Get cmake
# Get maven
# Get Python 2.7
# Get PCRE, which we need to build SWIG

# Get SWIG 3, which we need to parse C++11. This is going to be in Ubuntu Utopic
# and higher, but we have to work with the Ubuntu that the Java images use, so
# we have to build from source.
RUN wget http://downloads.sourceforge.net/swig/swig-3.0.5.tar.gz && \
    tar -xvzf swig-3.0.5.tar.gz && cd swig-3.0.5 && \
    ./configure --prefix=/usr && make && make install && cd .. && \
    rm -Rf swig-3.0.5 && rm swig-3.0.5.tar.gz

# Get SBT manually, since their HTTPS source segfaults my apt-get
RUN wget https://dl.bintray.com/sbt/debian/sbt-0.13.7.deb && \
    dpkg -i sbt-0.13.7.deb && rm sbt-0.13.7.deb

# Get the SDSL dependency and install it under /usr so we don't need to mess
# with compiler paths.
RUN git clone https://github.com/simongog/sdsl-lite.git && \
    cd sdsl-lite && ./install.sh /usr && cd .. && rm -Rf sdsl-lite

# Add in the entire Git repo that you are supposed to be running this from. We
# do this so we can build off a working copy, instead of going to Github and
# always building master.
ADD . sequence-graphs
# Make sure it is actually a git repo, make sure we have the submodules, and
# build.
RUN cd sequence-graphs && git init && git submodule init && \
    git submodule update && make
    
# Test our build
RUN cd sequence-graphs && make check

# Set everything up so we can enter the container sanely.
ENTRYPOINT ["/bin/bash"]
CMD ["-c", "cd sequence-graphs && /bin/bash"]
