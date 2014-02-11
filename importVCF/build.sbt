import AssemblyKeys._

assemblySettings

// Build a tool to import a VCF to sequence graph format.

packageArchetype.java_application

exportJars := true

name := "importVCF"

version := "0.1"

scalaVersion := "2.10.3"

// See unchecked warnings. See <http://stackoverflow.com/a/9417094/402891>
scalacOptions ++= Seq("-unchecked")

libraryDependencies += "org.apache.spark" %% "spark-core" % "0.9.0-incubating"

libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.9.0-incubating"

// We need command-line options

libraryDependencies += "org.rogach" %% "scallop" % "0.9.4"

// We get innovativemedicine's vcfimp VCF parser, which may not exist in any of
// our resolvers.
libraryDependencies += "ca.innovativemedicine" %% "vcfimp" % "0.7.0"

// We need a parallel VCF parser
libraryDependencies += "fi.tkk.ics.hadoop.bam" % "hadoop-bam" % "6.1-SNAPSHOT"

// We need log4j so we can tell Spark to be somewhat quieter and not swamp out
// output.
libraryDependencies += "log4j" % "log4j" % "1.2.17"

// We need some SLF4J logging implementation so it will stop logging complaints
// about not having a logging implementation. See
// <https://github.com/scalikejdbc/scalikejdbc/issues/21>.
// This implementation will pump log messages through log4j above.
libraryDependencies += "org.slf4j" % "slf4j-log4j12" % "1.7.5"

// We need an slf4j-api to go with our implementation
libraryDependencies += "org.slf4j" % "slf4j-api" % "1.7.5"


resolvers += "Akka Repository" at "http://repo.akka.io/releases/"

resolvers += "Hadoop-BAM" at "http://hadoop-bam.sourceforge.net/maven/"

resolvers += "Sonatype" at "http://oss.sonatype.org/content/repositories/snapshots/"

mainClass in assembly := Some("importVCF")

mergeStrategy in assembly := {
 case m if m.toLowerCase.endsWith("manifest.mf") => MergeStrategy.discard
 case m if m.toLowerCase.matches("meta-inf.*\\.sf$") => MergeStrategy.discard
 case "META-INF/services/org.apache.hadoop.fs.FileSystem" => MergeStrategy.concat
 case "reference.conf" => MergeStrategy.concat
 case "log4j.properties" => MergeStrategy.concat
 case _ => MergeStrategy.first
}

