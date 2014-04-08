import sbt._
import Keys._

object SGBuild extends Build {
    // This project holds the library and aggregates the subprojects for easy
    // building.
    lazy val sequencegraphs = Project(id = "sequencegraphs",
        base = file(".")) aggregate(
            LocalProject("importVCF"), 
            LocalProject("debug"), 
            LocalProject("exportVCF"),
            LocalProject("createIndex"),
            LocalProject("mapToIndex")
        )

    // This project depends on the library. Since it depends on a project that
    // aggregates it, we ned to use indirect lookup references with LocalProject
    // as described at <https://groups.google.com/forum/#!topic/simple-build-
    // tool/Yjc5gdejLvU>, instead of just passing the sequencegraphs and
    // importVCF objects to each other's definitions.
    lazy val importVCF = Project(id = "importVCF",
        base = file("importVCF")) dependsOn(
            // We depend on the root sequence graphs library
            LocalProject("sequencegraphs")
        )

    lazy val debug = Project(id = "debug",
        base = file("debug")) dependsOn(
            // We depend on the root sequence graphs library
            LocalProject("sequencegraphs")
        )
    lazy val exportVCF = Project(id = "exportVCF",
        base = file("exportVCF")) dependsOn(
            // We depend on the root sequence graphs library
            LocalProject("sequencegraphs")
        )
        
    lazy val createIndex = Project(id = "createIndex",
        base = file("createIndex")) dependsOn(
            // We depend on the root sequence graphs library
            LocalProject("sequencegraphs")
        )
        
    lazy val mapToIndex = Project(id = "mapToIndex",
        base = file("mapToIndex")) dependsOn(
            // We depend on the root sequence graphs library
            LocalProject("sequencegraphs")
        )

}
