##Dependencies

###To install Apache Spark with GraphX

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
./importVCF.sh <vcf file> <sample name> <output directory name>
```
