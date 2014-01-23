// Build a tool to dump a sequence graph as JSON

packageArchetype.java_application

exportJars := true

name := "debug"

version := "0.1"

scalaVersion := "2.9.3"

libraryDependencies += "org.apache.spark" %% "spark-core" % "0.8.1-incubating"

//libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.9.0-incubating"

// We need command-line options

libraryDependencies += "org.rogach" %% "scallop" % "0.9.4"

// We get innovativemedicine's vcfimp VCF parser, which may not exist in any of
// our resolvers. We need to specify a complicated version remapping function
// because vcfimp insists on building for Scala 2.9.2 only, which should be
// actually compatible with anything 2.9 but doesn't get picked up.
libraryDependencies += "ca.innovativemedicine" %% "vcfimp" % "0.7.0"

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"

resolvers += "Hadoop-BAM" at "http://hadoop-bam.sourceforge.net/maven/"

resolvers += "Sonatype" at "http://oss.sonatype.org/content/repositories/snapshots/"
