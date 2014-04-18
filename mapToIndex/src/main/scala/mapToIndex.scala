package edu.ucsc.genome.MapToIndex

import scala.collection.mutable.HashMap

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import org.apache.spark.rdd.RDD

import org.apache.spark.graphx
import org.apache.spark.graphx._

import edu.ucsc.genome._

import java.io.File
import java.nio.file.Paths

// We want to be able to loop over Java iterators: load a bunch of conversions.
// See <http://stackoverflow.com/a/1625252/402891>
import scala.collection.JavaConversions._

import scala.math._

// We want to parse command-line arguments
import org.rogach.scallop._

/**
 * Tool to create a ReferenceHierarchy from a set of FASTA files.
 */
object MapToIndex {
    def main(args: Array[String]) {
        
        // Option parsing code more or less stolen from various Scallop
        // examples.
        val opts = new ScallopConf(args) {
            guessOptionName = true
            banner("""Usage: mapToIndex [OPTION] index sequence
                |Map a sequence string to a reference hierarchy.
                |Options:
                |""".stripMargin)
            
            val cluster = opt[String](default = Some("local"),
                descr = "Run against this cluster URL")
            
            val index = trailArg[String](required = true,
                descr = "index path to save in")
                
            val pattern = trailArg[String](required = true,
                descr = "string to map")
            
            val version = opt[Boolean](noshort = true, 
                descr = "Print version")
            val help = opt[Boolean](noshort = true, 
                descr = "Show this message")

        } 
        
        // Set up logging
        ConsoleUtil.quietParquet
        
        // Set up Spark.
        
        // Set up serialization stuff for Spark so it can efficiently exchange
        // our Avro records.
        SequenceGraphKryoProperties.setupContextProperties()
        
        // Set the executors to use a more reasonable amount of memory.
        System.setProperty("spark.executor.memory", "25G")
        
        // Set the parallelism level to have enough reducers to not run out of
        // memory
        System.setProperty("spark.default.parallelism", "12")
        
        // The first thing we need is a Spark context. We would like to be able
        // to make one against any Spark URL: either "local" or soemthing like
        // "mesos://wherever.biz:1234". We need to feed it all the jars it needs
        // to run our code. Fortunately, they all live next to our jar, if the
        // sbt native packager has worked correctly.
        
        // What File is the jar that this class is from?
        val jarFile = new File(getClass.getProtectionDomain
            .getCodeSource.getLocation.toURI)
        
        // What files are in that directory (should all be .jars)? Make a list
        // of their string paths.
        val jarsToSend = jarFile.getParentFile.listFiles.map(_.toString).toSeq
            
        // Set up Spark, giving it the appropriate cluster URL, the SPARK_HOME
        // environment variable (which must be set!), and the list of jars we
        // have worked out.
        println("Initializing Spark")
        val sc = new SparkContext(opts.cluster.get.get, "createIndex", 
            System.getenv("SPARK_HOME"), jarsToSend)
        println("Spark initialized")
        
        // Get the hierarchy path
        val indexPath = opts.index.get.get
        
        // And the string to map
        val pattern = opts.pattern.get.get 

        // Make the ReferenceHierarchy from the index we built
        val hierarchy = new ReferenceHierarchy(sc, indexPath + "/hierarchy")
        
        println("Hierarchy loaded! Mapping...")
        
    }
}





















