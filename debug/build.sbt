import AssemblyKeys._

assemblySettings

packageArchetype.java_application

exportJars := true

name := "debug"

version := "0.1"

scalaVersion := "2.10.3"

// Don't generate path names that are too long for an encrypted home directory
// on Ubuntu. eCryptfs has a 143-character limit. See
// <https://bugs.launchpad.net/ubuntu/+source/linux/+bug/344878>
scalacOptions ++= Seq("-Xmax-classfile-name", "130")

libraryDependencies += "org.rogach" %% "scallop" % "0.9.4"

// We get innovativemedicine's vcfimp VCF parser, which may not exist in any of
// our resolvers.
libraryDependencies += "ca.innovativemedicine" %% "vcfimp" % "0.7.0"

libraryDependencies += "org.apache.spark" %% "spark-core" % "1.0.0-incubating-SNAPSHOT"

libraryDependencies += "org.apache.spark" %% "spark-graphx" % "1.0.0-incubating-SNAPSHOT"

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

// add hadoop bam and everything
libraryDependencies += "fi.tkk.ics.hadoop.bam" % "hadoop-bam" % "6.1"

libraryDependencies += "variant" % "variant" % "1.93"

libraryDependencies += "tribble" % "tribble" % "1.93"

libraryDependencies += "picard" % "picard" % "1.93"

libraryDependencies += "samtools" % "samtools" % "1.93"

libraryDependencies += "cofoja" % "cofoja" % "1.0"

mainClass in assembly := Some("debug")

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
