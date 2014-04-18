package edu.ucsc.genome

import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.graphx._

import org.apache.commons.io.FileUtils
import java.io._

import com.esotericsoftware.kryo._
import com.esotericsoftware.kryo.io.{Input, Output}
import com.twitter.chill.ScalaKryoInstantiator

import scala.reflect._

/**
 * Represents a reference hierarchy composed of reference structures at
 * different levels. Defined by a bottom-level index, a merging scheme for
 * building each level, and a possibly finished, or possibly null, graph of
 * Sides, Sites, Adjacencies, and Generalizations (with no negative vertex IDs).
 */
class ReferenceHierarchy(sc: SparkContext, var index: FMDIndex) {
    
    /**
     * Load a ReferenceHierarchy from the given path, using the given
     * SparkContext. The saved hierarchy must contain at least one level.
     */
    def this(sc: SparkContext, path: String) = {
        // Load an FMDIndex for a basename in that directory.
        this(sc, new FMDIndex(path + "/index.basename"))
    }
    
    
    /**
     * Save this ReferenceHierarchy to the specified directory.
     */
    def save(path: String) = {
        
        
    }
    
    
    
    
    
}


