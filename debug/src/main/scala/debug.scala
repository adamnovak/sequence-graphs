package edu.ucsc.genome.Debug

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import edu.ucsc.genome._

import java.io.File
import java.nio.file.Paths

// We want to read VCF files
import ca.innovativemedicine.vcf._
import ca.innovativemedicine.vcf.parsers._

// We want to be able to loop over Java iterators: load a bunch of conversions.
// See <http://stackoverflow.com/a/1625252/402891>
import scala.collection.JavaConversions._

// We want to parse command-line arguments
import org.rogach.scallop._;

import parquet.hadoop.ParquetInputFormat
import parquet.avro.{AvroParquetInputFormat, AvroReadSupport}
import parquet.hadoop.util.ContextUtil
import org.apache.hadoop.mapreduce.Job
import org.apache.spark.rdd.RDD

object Debug {
  def main(args: Array[String]) {
    
    // Option parsing code more or less stolen from various Scallop
    // examples.
    val opts = new ScallopConf(args) {
      guessOptionName = true
      banner("""Usage: debug [OPTION]... dirs
             |Dump sequence graph format in JSON.
             |Options:
             |""".stripMargin)
      
      val dirs = trailArg[String](required = true,
                                  descr = "Parquet/graph to import")
      
      val version = opt[Boolean](noshort = true, 
                                 descr = "Print version")
      val help = opt[Boolean](noshort = true, 
                              descr = "Show this message")
    }

    SequenceGraphKryoProperties.setupContextProperties()
    val sc = new SparkContext("local", "debug")
    val job = new Job(sc.hadoopConfiguration)

    val filePath = opts.dirs.get.get

    ParquetInputFormat.setReadSupportClass(job, classOf[AvroReadSupport[Side]])
    val sides: RDD[Side] = sc.newAPIHadoopFile(filePath + "/Sides",
                                               classOf[ParquetInputFormat[Side]], 
                                               classOf[Void], classOf[Side],
                                               ContextUtil.getConfiguration(job)).map(p => p._2)
    println (sides.count() + " Sides:")
    sides.foreach(println(_))
    
    ParquetInputFormat.setReadSupportClass(job, classOf[AvroReadSupport[AlleleGroup]])
    val ag = sc.newAPIHadoopFile(filePath + "/AlleleGroups",
                                 classOf[ParquetInputFormat[AlleleGroup]], 
                                 classOf[Void], classOf[AlleleGroup],
                                 ContextUtil.getConfiguration(job)).map(p => p._2)
    println (ag.count() + " Allele Groups:")
    ag.foreach(println(_))
    
    ParquetInputFormat.setReadSupportClass(job, classOf[AvroReadSupport[Adjacency]])
    val ad = sc.newAPIHadoopFile(filePath + "/Adjacencies",
                                 classOf[ParquetInputFormat[Adjacency]], 
                                 classOf[Void], classOf[Adjacency],
                                 ContextUtil.getConfiguration(job)).map(p => p._2)
    println (ad.count() + " Adjacencies:")
    ad.foreach(println(_))

    ParquetInputFormat.setReadSupportClass(job, classOf[AvroReadSupport[Anchor]])
    val an = sc.newAPIHadoopFile(filePath + "/Anchors",
                                 classOf[ParquetInputFormat[Anchor]], 
                                 classOf[Void], classOf[Anchor],
                                 ContextUtil.getConfiguration(job)).map(p => p._2)
    println (an.count() + " Anchors:")
    an.foreach(println(_))
  }
}





















