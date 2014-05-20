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
                descr = "cluster URL to run against")
            
            val index = trailArg[String](required = true,
                descr = "index path to load")
                
            val pattern = trailArg[String](required = true, 
                descr = "string to map")
            
            val repeat = opt[Int](noshort = true, default=Some(1),
                descr = "number of times to repeat the mapping")
            
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
        
        // Get the hierarchy path
        val indexPath = opts.index.get.get
        
        // And the string to map
        val pattern = opts.pattern.get.get 

        // Make the ReferenceHierarchy from the index we built
        val hierarchy = new ReferenceHierarchy(indexPath)
        
        println("Hierarchy loaded! Mapping %s times..."
            .format(opts.repeat.get.get))
        
        // We want to track time per mapping. See
        // <http://rosettacode.org/wiki/Time_a_function#Scala>
        val startTime = System.currentTimeMillis
        
        // Map to all levels
        var mappings: Seq[Seq[Option[Side]]] = null
        
        for(i <- 0 until opts.repeat.get.get) {
            // Repeat the mapping several times for benchmarking purposes.
            mappings = hierarchy.map(pattern)
        }
        
        // Stop the timer
        val endTime = System.currentTimeMillis
        
        // Get the left and right mappings
        val leftMappings = hierarchy.map(pattern, Face.LEFT)
        
        val rightMappings = hierarchy.map(pattern, Face.RIGHT)
        
        // Have a little function to format the sides
        def mappingToString(m: Option[Side]): String = {
            m match {
                case Some(side) => 
                    "%d%s".format(side.coordinate, side.face match {
                        case Face.LEFT => "L"
                        case Face.RIGHT => "R"
                    })
                case None => "          "
            }
        }
        
        // Bind up all the mappings together
        val allMappings = (leftMappings, mappings, rightMappings).zipped
        
        for((mappingTuple, index) <- allMappings.toSeq.zipWithIndex) {
            println("Mappings to level %d:".format(index))
            
            println(Seq("Base", "Left", "Both", "Right").map(_.padTo(10, ' '))
                .mkString("\t"))
            
            // Get the individual mapping lists
            val (left, center, right) = mappingTuple
            
            // Zip up with the pattern. Not really a zip; we use "transpose" to
            // flip which direction is the outer list.
            val rows = List(pattern.toSeq.map(_.toString),
                left.map(mappingToString _), center.map(mappingToString _), 
                right.map(mappingToString _))
                .transpose
            
            rows.map { row =>
                // Print each row
                println(row.map(_.padTo(10, ' ')).mkString("\t"))
            }
        }
        
        
        // Report time per mapping.
        println("map() took %f ms/call".format((endTime - startTime).toDouble / 
            opts.repeat.get.get))
        
    }
}





















