package edu.ucsc.genome.ExportVCF

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import edu.ucsc.genome.{SequenceGraph, AlleleGroup, Side, Adjacency, Anchor, Edge, Allele=>SGAllele, Position, SequenceGraphKryoProperties}
import org.broadinstitute.variant.variantcontext.{Allele, Genotype, VariantContext, GenotypeBuilder, GenotypesContext, VariantContextBuilder}
import fi.tkk.ics.hadoop.bam.VariantContextWritable
import parquet.hadoop.ParquetInputFormat
import parquet.avro.{AvroParquetInputFormat, AvroReadSupport}
import parquet.hadoop.util.ContextUtil
import org.apache.hadoop.mapreduce.Job
import org.apache.spark.rdd.RDD
import org.apache.hadoop.io.LongWritable
import scala.collection.JavaConversions._

import java.io.File
import java.nio.file.Paths


// We want to be able to loop over Java iterators: load a bunch of conversions.
// See <http://stackoverflow.com/a/1625252/402891>
import scala.collection.JavaConversions._

// We want to parse command-line arguments
import org.rogach.scallop._;

object ExportVCF {
  def main(args: Array[String]) {
    
    // Option parsing code more or less stolen from various Scallop
    // examples.
    val opts = new ScallopConf(args) {
      guessOptionName = true
      banner("""Usage: exportVCF [OPTION]... inputDir vcfFile sample
             |Export a sequence graph to VCF file.
             |Options:
             |""".stripMargin)
             
      val cluster = opt[String](default = Some("local"),
                descr = "Run against this cluster URL") 
      
      val directory = trailArg[String](required = true,
                                       descr = "Directory to read from")
      val vcfFile = trailArg[String](required = true,
                                     descr = "VCF file to save in")
      val sample = trailArg[String](required = true,
                                    descr = "Sample name")
      val version = opt[Boolean](noshort = true, 
                                 descr = "Print version")
      val help = opt[Boolean](noshort = true, 
                              descr = "Show this message")
    } 

    val export = new ExportVCF(opts.cluster.get.get, opts.directory.get.get, 
        opts.vcfFile.get.get, opts.sample.get.get)

    export.export()
  }
}

class ExportVCF (cluster: String, directory: String, vcfFile: String, 
    sample: String) extends Serializable {

  def export() {

    // Set up serialization stuff for Spark so it can efficiently exchange our
    // Avro records.
    SequenceGraphKryoProperties.setupContextProperties()

    // The first thing we need is a Spark context. We would like to be able to
    // make one against any Spark URL: either "local" or soemthing like
    // "mesos://wherever.biz:1234".
    
    // Unfortunately, when we use a Mesos URL, Spark dies. In an effor to
    // appease it, we try to give it all of the .jars our program uses.
    // Unfortunately, in general a Java program has no idea where to find these
    // jars. However, this program is built with the sbt native packager, which
    // puts all the jars retrieved from various repositories and from subproject
    // dependencies together in a big lib folder. So we can look to see what jar
    // this class has been loaded from, and load up all the jars in that
    // directory.
    
    // What File is the jar that this class is from?
    val jarFile = new File(classOf[ExportVCF].getProtectionDomain.getCodeSource
        .getLocation.toURI)
    
    // What files are in that directory (should all be .jars)? Make a list of
    // their string paths.
    val jarsToSend = jarFile.getParentFile.listFiles.map(_.toString).toSeq
        
    // Set up Spark, giving it the appropriate cluster URL, the SPARK_HOME
    // environment variable if set, and the list of jars we have worked out.
    println("Initializing Spark")
    val sc = new SparkContext(cluster, "exportVCF", 
        System.getenv("SPARK_HOME"), jarsToSend)
    println("Spark initialized")
    val job = new Job(sc.hadoopConfiguration)

    // read data in
    val filePath = directory


    // Load sequence graph Sides.
    val sides: RDD[Side] = SequenceGraph.readRDDFromParquet(sc,
        filePath + "/Sides")
        
    // And AlleleGroups
    val ag: RDD[AlleleGroup] = SequenceGraph.readRDDFromParquet(sc,
        filePath + "/AlleleGroups")

    // And Adjacencies
    val ad: RDD[Adjacency] = SequenceGraph.readRDDFromParquet(sc,
        filePath + "/Adjacencies")

    // And Anchors
    val an: RDD[Anchor] = SequenceGraph.readRDDFromParquet(sc,
        filePath + "/Anchors")
 
    // get phasing sets
    // key all sides by their ids
    val keyedSides = sides.keyBy(_.getId).cache()
    // for each start anchor, key with the id of the left side, and join with the side
    val startAnch = an.keyBy(a => a.getEdge.getLeft)
      .join(keyedSides)
      .map(kv => kv._2._2.getPosition)
    // for each end anchor, key with the id of the left side, and join with the side
    val endAnch = an.keyBy(a => a.getEdge.getRight)
      .join(keyedSides)
      .map(kv => kv._2._2.getPosition)
    // to recover phasing sets, zip all anchors to discover regions that are phased, then
    // get counts.
    // regions that are phased will have multiple anchors show up at a single location
    // and will have counts greater than 1
    // we then group all anchors by contig, and build a phase set
    val phaseSets = startAnch.zip(endAnch)
      .map(r => (r, 1))
      .reduceByKey(_ + _)
      .groupBy(_._1._1.getContig)
      .flatMap(kv => {
        val (contig, ranges) = kv
        val builder = new PhaseSetBuilder(contig)
        
        ranges.map(t => {
          val ((s, e), c) = t

          // flip tuple around and get base positions on contig
          (c, s.getBase, e.getBase)
        }).map(kv => (kv._1, kv._2.toLong, kv._3.toLong))
        .foreach(builder.add(_)) // loop over and add anchors to phase set

        // build phase sets for this contig
        builder.toSets()
      })
      .collect()

    // set contig and sample info in header
    SequenceGraphVCFOutputFormat.addSample(sample)
    
    // add contigs to vcf header:
    // map over sides to get position value, then get contig and turn into string and uniquify.
    // then, pass to picard via hadoop-bam output format wrapper
    sides.map(_.getPosition)
      .map(_.getContig)
      .map(_.toString)
      .distinct()
      .map(new java.lang.String(_))
      .collect()
      .foreach(SequenceGraphVCFOutputFormat.addContig(_))
    
    // get allele info:
    // join with sides to recover position values
    val alleles = ag.keyBy(a => a.getEdge.getLeft)
      .join(keyedSides)
      .map(kv => (kv._2._2.getPosition, kv._2._1))
      .groupByKey()
      
    // convert to variants
    val variantContexts: RDD[(LongWritable,VariantContextWritable)] = alleles.map(kv => convertToVariantContext(kv._1, 
                                                                                                                kv._2.toList, 
                                                                                                                phaseSets, 
                                                                                                                sample))
      .keyBy(vcw => new LongWritable(vcw.get.getStart.toLong))

    // write out
    variantContexts.saveAsNewAPIHadoopFile(vcfFile,
                                           classOf[LongWritable], 
                                           classOf[VariantContextWritable],
                                           classOf[SequenceGraphVCFOutputFormat[LongWritable]])

    // TODO: Flush logging output from Parquet/Spark stuff.
    println("VCF exported")
  }

  /**
   * Converts data to a Picard variant context.
   *
   * @param position Reference position of this polymorphism.
   * @param alleles All alleles which segregate at this site
   * @param phaseSets A list of all phase sets in this callset.
   * @param sample Sample name.
   * @return A Picard variant context.
   */
  def convertToVariantContext(position: Position,
                              alleles: List[AlleleGroup],
                              phaseSets: Array[PhaseSet],
                              sample: String): VariantContextWritable = {
    // convert alleles to Picard allele
    val alleleNoGroup = alleles.flatMap((r: AlleleGroup) => r.getAllele() match {
      case (a: SGAllele) => { // we have an allele
        var l = List[SGAllele]()

        // get ploidy lower bound
        val ploidy = r.getPloidy.getLower()
        
        // create ploidy count of alleles
        for (i <- 0 until ploidy) {
          l = a :: l
        }

        l
      }
      case _ => List[SGAllele]()
    })

    // get reference at site of variation
    val refBases = alleleNoGroup.head.getReferenceBases

    // turn alleles into a set of bases
    val alleleBases = alleleNoGroup.map(_.getBases())

    /**
     * Turns a string of bases into a Picard Allele. Sets whether allele is the
     * reference allele.
     *
     * @param bases Allele as a string of bases.
     * @return A Picard allele corresponding to this string of bases.
     */
    def buildAllele (bases: String): Allele = {
      val ref = (refBases == bases)
      
      Allele.create(bases, ref)
    }

    // alleles fed to variant context need to be unique, genotypes do not have this constraint
    val alleleConvWoRef = alleleBases.map(buildAllele)
    val alleleConvVC = alleleBases.distinct.map(buildAllele)

    // add reference allele if it doesn't exist in set
    val alleleConv = if(alleleConvVC.map(_.isReference).reduce(_ || _)) {
      alleleConvVC
    } else {
      buildAllele(refBases) :: alleleConvVC
    }
    
    // build genotypes
    val gb = new GenotypeBuilder(sample, alleleConvWoRef)
    
    // check for phasing:
    // if phase set contains a phase set who is:
    // - on the same contig
    // - starts before and ends after
    // our current position, then we are phased
    val contig = position.getContig()
    val pos = position.getBase()
    val matchingPhaseSet = phaseSets.filter(v => v.contig == contig)
          .filter(v => v.start < pos)
          .filter(v => v.end > pos)
    if (matchingPhaseSet.length == 1) {
      gb.phased(true)
      gb.attribute("PS", matchingPhaseSet.head.start)
    }
    
    // finish making genotypes
    val genotype = gb.make()

    // make genotype context... - ugly manipulation...
    val genotypeList: java.util.List[Genotype] = bufferAsJavaList(List(genotype).toBuffer)
    val genotypeArrayList: java.util.ArrayList[Genotype] = new java.util.ArrayList(genotypeList)
    val sampleMap: java.util.Map[java.lang.String, java.lang.Integer] = mapAsJavaMap(Map(new java.lang.String(sample) -> new java.lang.Integer(0)))
    val sampleOrder: java.util.List[String] = bufferAsJavaList(List(new java.lang.String(sample)).toBuffer)
    val gc = GenotypesContext.create(genotypeArrayList, sampleMap, sampleOrder)

    // build variant context
    val vcb = new VariantContextBuilder("SequenceGraphExport",
                                        position.contig,
                                        position.base,
                                        position.base + refBases.length - 1,
                                        alleleConv)
    vcb.genotypes(gc)
    val vc = vcb.make()

    // make and return writable variant context
    val vcw = new VariantContextWritable
    vcw.set(vc)

    vcw
  }
}





















