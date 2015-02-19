# This file can be used by Docker to create a container within which the
# sequence graph test code can run successfully. It can also double as
# installation instructions.

# Get Ubuntu Whatever with a JDK on it, and up to date packages.
FROM dockerfile/java:oracle-java7
MAINTAINER Adam Novak <anovak@soe.ucsc.edu>
RUN apt-get update
RUN apt-get upgrade -y

# Get Git
RUN apt-get install -y git

# Put in GCC. Assume we get 4.8+, since we can't force a version on the java
# base image.
RUN apt-get install -y build-essential gcc

# Get Boost
RUN apt-get install -y libboost-all-dev

# Get sparse hash library
RUN apt-get install -y libsparsehash-dev

# Get cppunit
RUN apt-get install -y libcppunit-dev

# Get Swig
RUN apt-get install -y swig

# Get cmake
RUN apt-get install -y cmake

# Get SBT manually, since their HTTPS source segfaults my apt-get
RUN wget https://dl.bintray.com/sbt/debian/sbt-0.13.7.deb
RUN dpkg -i sbt-0.13.7.deb

# Get maven
RUN apt-get install -y maven

# Get Python 2.7
RUN apt-get install -y python2.7

# Get the SDSL dependency and install it under /usr so we don't need to mess
# with compiler paths.
RUN git clone https://github.com/simongog/sdsl-lite.git && cd sdsl-lite && ./install.sh /usr

# Clone my code and all the submodules
RUN git clone --recursive https://github.com/adamnovak/sequence-graphs.git && cd sequence-graphs && make

# Set everything up so we can use the container like a command
ENTRYPOINT ["/sequence-graphs/createIndex/createIndex"]
CMD ["--help"]
