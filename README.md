##Dependencies

###To install Spark 0.9.0 with Mesos support

Support for the Mesos cluster manager in anything other than coarse scheduling mode is broken in Spark 0.9.0. In order to use Mesos in a reasonable scheduling mode, it is necessary to install a version of Spark with a patch that has not yet been merged into the main development repository.

Make sure to set `SPARK_HADOOP_VERSION=<your Hadoop version>` and `SPARK_YARN=<true or false depending on whether you want to use the new map-reduce>`.

```
export SPARK_HADOOP_VERSION=<your Hadoop version>
export SPARK_YARN=<whether you want to use YARN or not>

git clone git@github.com:apache/incubator-spark.git
cd incubator-spark
git fetch git@github.com:bijaybisht/incubator-spark.git SPARK-1052
git checkout -b SPARK-1052 FETCH_HEAD
sbt clean
sbt assembly
sbt publish-local
```

### To install ADAM

```
git clone git@github.com:bigdatagenomics/adam.git
cd adam
```

Make sure to edit `pom.xml`, following the instructions to fill in your Hadoop version. If you are running against a Mesos cluster, it is vitally important to be correct about whether Hadoop 2.2.0 or greater is being used or not, because this determines the version of protobuf that gets used. An incorrect version of protobuf can lead to jobs segfaulting when trying to load Mesos

```
MAVEN_OPTS="-Xmx512m -XX:MaxPermSize=128m" mvn clean package install
```

##Installation

###Building

```
sbt stage
```

###Testing

```
sbt test
```

###Packaging

A single, JAR which includes all dependencies can be created by running:

```
sbt assembly
```

###Running command-line tools

```
./importVCF.sh --parquet-dir <output directory absolute path> <vcf file> <sample name>
./importVCF.sh --cluster mesos://localhost:5050 --parquet-dir <output directory absolute path> <vcf file> <sample name> 
./importVCF.sh --dot-file <output file name> <vcf file> <sample name>
```
