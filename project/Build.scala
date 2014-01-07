import sbt._
import Keys._

object SGBuild extends Build {
    // This project holds the library and aggregates the subprojects for easy
    // building.
    lazy val sequencegraphs = Project(id = "sequencegraphs",
        base = file(".")) aggregate(LocalProject("importVCF"))

    // This project depends on the library. Since it depends on a project that
    // aggregates it, we ned to use indirect lookup references with LocalProject
    // as described at <https://groups.google.com/forum/#!topic/simple-build-
    // tool/Yjc5gdejLvU>, instead of just passing the sequencegraphs and
    // importVCF objects to each other's definitions.
    lazy val importVCF = Project(id = "importVCF",
        base = file("importVCF")) dependsOn(LocalProject("sequencegraphs"))

}
