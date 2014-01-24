// Build a tool to dump a sequence graph as VCF

packageArchetype.java_application

exportJars := true

name := "exportVCF"

version := "0.1"

scalaVersion := "2.9.3"

libraryDependencies += "org.apache.spark" %% "spark-core" % "0.8.1-incubating"

//libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.9.0-incubating"

// We need command-line options

libraryDependencies += "org.rogach" %% "scallop" % "0.9.4"

// add hadoop bam and everything

libraryDependencies += "fi.tkk.ics.hadoop.bam" % "hadoop-bam" % "6.1-SNAPSHOT"

libraryDependencies += "variant" % "variant" % "1.93"

libraryDependencies += "tribble" % "tribble" % "1.93"

libraryDependencies += "picard" % "picard" % "1.93"

libraryDependencies += "samtools" % "samtools" % "1.93"

libraryDependencies += "cofoja" % "cofoja" % "1.0"

// Get mesos ourselves

libraryDependencies += "org.apache.mesos" % "mesos" % "0.13.0"

// add adam

libraryDependencies += "edu.berkeley.cs.amplab.adam" % "adam-core" % "0.6.1-SNAPSHOT"

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"

resolvers += "Hadoop-BAM" at "http://hadoop-bam.sourceforge.net/maven/"

resolvers += "Sonatype" at "http://oss.sonatype.org/content/repositories/snapshots/"
