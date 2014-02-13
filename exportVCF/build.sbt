
import AssemblyKeys._

assemblySettings

// Build a tool to import a VCF to sequence graph format.

packageArchetype.java_application

exportJars := true

name := "exportVCF"

scalaVersion := "2.10.3"

// Don't generate path names that are too long for an encrypted home directory
// on Ubuntu. eCryptfs has a 143-character limit. See
// <https://bugs.launchpad.net/ubuntu/+source/linux/+bug/344878>
scalacOptions ++= Seq("-Xmax-classfile-name", "130")

libraryDependencies += "org.apache.spark" %% "spark-core" % "0.9.0-incubating"

libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.9.0-incubating"

// We need this library for testing

libraryDependencies += "org.scalatest" %% "scalatest" % "1.9.2" % "test"

// We need to be able to write out graph pictures

libraryDependencies += "org.kohsuke" % "graphviz-api" % "1.1"

// We need to be able to write Avro in Parquet

libraryDependencies += "com.twitter" % "parquet-avro" % "1.3.2"

libraryDependencies += "org.rogach" %% "scallop" % "0.9.4"

// We need slf4j for logging

libraryDependencies += "org.slf4j" % "slf4j-api" % "1.7.5"

// We need Apache Commons Lang because it's useful (for escaping strings in .dot
// files)
libraryDependencies += "org.apache.commons" % "commons-lang3" % "3.2.1"

// We need Apache Commons IO because "rm -Rf" takes about a hunderd lines of
// Java.
libraryDependencies += "commons-io" % "commons-io" % "2.4"

// add hadoop bam and everything
libraryDependencies += "fi.tkk.ics.hadoop.bam" % "hadoop-bam" % "6.1-SNAPSHOT"

libraryDependencies += "variant" % "variant" % "1.93"

libraryDependencies += "tribble" % "tribble" % "1.93"

libraryDependencies += "picard" % "picard" % "1.93"

libraryDependencies += "samtools" % "samtools" % "1.93"

libraryDependencies += "cofoja" % "cofoja" % "1.0"

// Add the local Maven repository as a resolver, so we can use development ADAM
// versions installed with "mvn install". This is in turn necessary in order to
// use any Hadoop version other than the one that the official ADAM artifacts
// are built against. See <http://stackoverflow.com/a/10778151/402891>

resolvers += "Local Maven" at Path.userHome.asFile.toURI.toURL + ".m2/repository"

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"

resolvers += "Hadoop-BAM" at "http://hadoop-bam.sourceforge.net/maven/"

resolvers += "Sonatype" at "http://oss.sonatype.org/content/repositories/snapshots/"

// Set externalResolvers up ourselves, so that the default Maven Central
// repository becomes the last resolver instead of the second (after Ivy cache).
// This allows us to override things in Maven Central with those in our local
// Maven repository, which we need to do in order to build against a version of
// ADAM that itself has been built against the correct version of Hadoop.

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"

mainClass in assembly := Some("exportVCF")

mergeStrategy in assembly := {
 case m if m.toLowerCase.endsWith("manifest.mf") => MergeStrategy.discard
 case m if m.toLowerCase.matches("meta-inf.*\\.sf$") => MergeStrategy.discard
 case "META-INF/services/org.apache.hadoop.fs.FileSystem" => MergeStrategy.concat
 case "reference.conf" => MergeStrategy.concat
 case "log4j.properties" => MergeStrategy.concat
 case _ => MergeStrategy.first
}
