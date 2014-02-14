package edu.ucsc.genome

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

/**
 * Companion object for the SparkSuite trait, containing the "static" things
 * that need to be set up only once.
 */
object SparkSuite {
    var sparkContext: Option[SparkContext] = None
    
    def getSpark: SparkContext = {
        this.synchronized {
            // Make sure we seriously only ever do this stuff once. Synchronize
            // on a singleton.
            sparkContext match { 
                case None =>
                    // Do our static initialization
                    
                    // Make sure Kryo is on.
                    SequenceGraphKryoProperties.setupContextProperties()
                    
                    // Set the SparkContext
                    sparkContext = Some(new SparkContext("local", "SparkTest"))
                    // Recurse
                    getSpark
                // If we already made a SparkContext, use that.
                case Some(context) => context
            }
        }
    }
}

/**
 * A trait for tests that need to use Spark. Makes sure Spark is set up only
 * once, and makes sure it is in scope.
 */
trait SparkSuite {
    import SparkSuite._
    
    // Grab the SparkContext from the companion object and pop it into scope.
    // Also does any static initialization needed to make Spark work right.
    val sc = getSpark

}
