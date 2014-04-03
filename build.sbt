import AssemblyKeys._

assemblySettings

// Project info
name := "Sequence Graph API"

version := "0.1"

// Set up the Avro code generator
seq( sbtavro.SbtAvro.avroSettings : _*)

(stringType in avroConfig) := "String"

(version in avroConfig) := "1.7.5"

// Set up dependency graph drawing
net.virtualvoid.sbt.graph.Plugin.graphSettings

// This is a library, not an application, so no native packager setup is needed.

scalaVersion := "2.10.3"

// Don't generate path names that are too long for an encrypted home directory
// on Ubuntu. eCryptfs has a 143-character limit. See
// <https://bugs.launchpad.net/ubuntu/+source/linux/+bug/344878>
scalacOptions ++= Seq("-Xmax-classfile-name", "130")

// We need Spark and GraphX
libraryDependencies += "org.apache.spark" %% "spark-core" % "0.9.0-incubating"

libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.9.0-incubating"

// We need RLCSA (with appropriate natives) for mapping. Install from
// https://github.com/adamnovak/rlcsa
libraryDependencies += "fi.helsinki.cs" % "rlcsa" % "1.0.0-SNAPSHOT"

// We need this library for testing
libraryDependencies += "org.scalatest" %% "scalatest" % "1.9.2" % "test"

// We need to be able to write out graph pictures
libraryDependencies += "org.kohsuke" % "graphviz-api" % "1.1"

// We need to be able to write Avro in Parquet
libraryDependencies += "com.twitter" % "parquet-avro" % "1.3.2"

// We need slf4j for logging
libraryDependencies += "org.slf4j" % "slf4j-api" % "1.7.5"

// We need Apache Commons Lang because it's useful (for escaping strings in .dot
// files)
libraryDependencies += "org.apache.commons" % "commons-lang3" % "3.2.1"

// We need Apache Commons IO because "rm -Rf" takes about a hunderd lines of
// Java.
libraryDependencies += "commons-io" % "commons-io" % "2.4"

// We need Google Guava for its UnionFind disjoint sets data structure.
libraryDependencies += "com.google.guava" % "guava" % "16.0.1"

// add hadoop bam and everything
libraryDependencies += "fi.tkk.ics.hadoop.bam" % "hadoop-bam" % "6.1"

// We need FASTA reading capability.
libraryDependencies += "org.biojava" % "biojava3-core" % "3.0.6"

libraryDependencies += "variant" % "variant" % "1.93"

libraryDependencies += "tribble" % "tribble" % "1.93"

libraryDependencies += "picard" % "picard" % "1.93"

libraryDependencies += "samtools" % "samtools" % "1.93"

libraryDependencies += "cofoja" % "cofoja" % "1.0"

mergeStrategy in assembly := {
 case m if m.toLowerCase.endsWith("manifest.mf") => MergeStrategy.discard
 case m if m.toLowerCase.matches("meta-inf.*\\.sf$") => MergeStrategy.discard
 case "META-INF/services/org.apache.hadoop.fs.FileSystem" => MergeStrategy.concat
 case "reference.conf" => MergeStrategy.concat
 case "log4j.properties" => MergeStrategy.concat
 case _ => MergeStrategy.first
}

// RLCSA comes from the local Maven repo. See
// <http://stackoverflow.com/a/10778151/402891>
resolvers += "Local Maven" at Path.userHome.asFile.toURI.toURL + ".m2/repository"

resolvers += "Sonatype" at "http://oss.sonatype.org/content/repositories/snapshots/"

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"

resolvers += "Hadoop-BAM" at "http://hadoop-bam.sourceforge.net/maven/"

resolvers += "BioJava repository" at "http://www.biojava.org/download/maven/"
