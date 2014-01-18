// Build a tool to import a VCF to sequence graph format.

packageArchetype.java_application

exportJars := true

name := "importVCF"

version := "0.1"

scalaVersion := "2.10.3"

libraryDependencies += "org.apache.spark" %% "spark-core" % "0.9.0-incubating"

libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.9.0-incubating"

// We need command-line options

libraryDependencies += "org.rogach" %% "scallop" % "0.9.4"

// We get innovativemedicine's vcfimp VCF parser, which may not exist in any of
// our resolvers.
libraryDependencies += "ca.innovativemedicine" %% "vcfimp" % "0.7.0"

// We need log4j so we can tell Spark to be somewhat quieter and not swamp out
// output.
libraryDependencies += "log4j" % "log4j" % "1.2.17"

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"
