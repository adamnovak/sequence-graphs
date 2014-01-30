##Dependencies

###To install Apache Spark

If you are going to use Spark to run executors for a cluster, make sure to set `SPARK_HADOOP_VERSION=<your Hadoop version>` and `SPARK_YARN=<true or false depending on whether you want to use the new map-reduce>` first.

```
git clone https://github.com/apache/incubator-spark
cd incubator-spark
git checkout v0.8.1-incubating
sbt publish-local
```

###To install vcfimp

```
git clone https://github.com/adamnovak/vcfimp.git
cd vcfimp
sbt "+ vcfimp/publish-local"
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

###Running command-line tools

```
./importVCF.sh --parquet-dir <output directory name> <vcf file> <sample name> 
./importVCF.sh --dot-file <output file name> <vcf file> <sample name>
```
