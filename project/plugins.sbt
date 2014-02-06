// Add the sbt Avro code generation plugin.

addSbtPlugin("com.cavorite" % "sbt-avro" % "0.3.2")

// Add the native packager for sbt-free launch scripts.

addSbtPlugin("com.typesafe.sbt" % "sbt-native-packager" % "0.6.4")

// Add a cool plugin so we can debug transitive dependencies visually
addSbtPlugin("net.virtual-void" % "sbt-dependency-graph" % "0.7.4")
