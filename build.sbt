// Set up the Avro code generator

seq( sbtavro.SbtAvro.avroSettings : _*)

(stringType in avroConfig) := "String"

// This is a library, not an application, so no native packager setup is needed.

name := "Sequence Graph API"

version := "0.1"

scalaVersion := "2.9.3"

libraryDependencies += "org.apache.spark" %% "spark-core" % "0.8.1-incubating"

//libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.8.1-incubating"

// We need this library for testing

libraryDependencies += "org.scalatest" %% "scalatest" % "1.9.2" % "test"

// We need to be able to write out graph pictures

libraryDependencies += "org.kohsuke" % "graphviz-api" % "1.1"

// We need to be able to write Avro in Parquet

libraryDependencies += "com.twitter" % "parquet-avro" % "1.3.2"

// We need slf4j for logging

libraryDependencies += "org.slf4j" % "slf4j-api" % "1.7.2"

// We need Apaceh Commons Lang because it's useful (for escaping strings in .dot
// files)
libraryDependencies += "org.apache.commons" % "commons-lang3" % "3.2.1"

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"

