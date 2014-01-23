##Dependencies

###To install Apache Spark

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
./importVCF.sh --parquet-dir <output directory name> <vcf file> <sample name> 
./importVCF.sh --dot-file <output file name> <vcf file> <sample name>
```
