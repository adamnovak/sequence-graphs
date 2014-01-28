name := "Simple Project"

version := "1.0"

scalaVersion := "2.9.3"

packageArchetype.java_application

libraryDependencies += "org.apache.spark" %% "spark-core" % "0.8.1-incubating"

//libraryDependencies += "org.apache.spark" %% "spark-core" % "0.9.0-incubating-SNAPSHOT"

//libraryDependencies += "org.apache.spark" %% "spark-graphx" % "0.9.0-incubating-SNAPSHOT"

//libraryDependencies += "org.apache.mesos" % "mesos" % "0.13.0"

resolvers += "Akka Repository" at "http://repo.akka.io/releases/"

resolvers += "Hadoop-BAM" at "http://hadoop-bam.sourceforge.net/maven/"

resolvers += "Sonatype" at "http://oss.sonatype.org/content/repositories/snapshots/"

