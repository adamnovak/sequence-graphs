// Build a tool to import a VCF to sequence graph format.

packageArchetype.java_application

exportJars := true

name := "importVCF"

version := "0.1"

scalaVersion := "2.9.3"

libraryDependencies += "org.apache.spark" %% "spark-core" % "0.9.0-incubating-SNAPSHOT"

libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.9.0-incubating-SNAPSHOT"

// We need command-line options

libraryDependencies += "org.rogach" %% "scallop" % "0.9.4"

// We get innovativemedicine's vcfimp VCF parser, which may not exist in any of
// our resolvers. We need to specify a complicated version remapping function
// because vcfimp insists on building for Scala 2.9.2 only, which should be
// actually compatible with anything 2.9 but doesn't get picked up.
libraryDependencies += "ca.innovativemedicine" % "vcfimp" % "0.6.1" cross 
    CrossVersion.binaryMapped {
        case "2.9.3" => "2.9.2" 
        case "2.9.1" => "2.9.2" 
        case "2.9.0" => "2.9.2"
        case x => x
    }

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"
