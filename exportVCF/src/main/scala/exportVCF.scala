package edu.ucsc.genome.ExportVCF

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import edu.ucsc.genome.{AlleleGroup, Side, Adjacency, Anchor, Edge, Allele=>SGAllele, Position, SequenceGraphKryoProperties}
import org.broadinstitute.variant.variantcontext.{Allele, Genotype, VariantContext, GenotypeBuilder, GenotypesContext, VariantContextBuilder}
import fi.tkk.ics.hadoop.bam.VariantContextWritable
import parquet.hadoop.ParquetInputFormat
import parquet.avro.{AvroParquetInputFormat, AvroReadSupport}
import parquet.hadoop.util.ContextUtil
import org.apache.hadoop.mapreduce.Job
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._
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

    val export = new ExportVCF(opts.directory.get.get, opts.vcfFile.get.get, opts.sample.get.get)

    export.export()
  }
}

class ExportVCF (directory: String, vcfFile: String, sample: String) extends Serializable {

  def export() {

    // read data in
    SequenceGraphKryoProperties.setupContextProperties()
    val sc = new SparkContext("local", "exportVCF")
    val job = new Job(sc.hadoopConfiguration)

    val filePath = directory

    ParquetInputFormat.setReadSupportClass(job, classOf[AvroReadSupport[Side]])
    val sides: RDD[Side] = sc.newAPIHadoopFile(filePath + "/Sides",
                                               classOf[ParquetInputFormat[Side]], 
                                               classOf[Void], classOf[Side],
                                               ContextUtil.getConfiguration(job)).map(p => p._2)
    
    ParquetInputFormat.setReadSupportClass(job, classOf[AvroReadSupport[AlleleGroup]])
    val ag = sc.newAPIHadoopFile(filePath + "/AlleleGroups",
                                 classOf[ParquetInputFormat[AlleleGroup]], 
                                 classOf[Void], classOf[AlleleGroup],
                                 ContextUtil.getConfiguration(job)).map(p => p._2)
    
    ParquetInputFormat.setReadSupportClass(job, classOf[AvroReadSupport[Adjacency]])
    val ad = sc.newAPIHadoopFile(filePath + "/Adjacencies",
                                 classOf[ParquetInputFormat[Adjacency]], 
                                 classOf[Void], classOf[Adjacency],
                                 ContextUtil.getConfiguration(job)).map(p => p._2)

    ParquetInputFormat.setReadSupportClass(job, classOf[AvroReadSupport[Anchor]])
    val an = sc.newAPIHadoopFile(filePath + "/Anchors",
                                 classOf[ParquetInputFormat[Anchor]], 
                                 classOf[Void], classOf[Anchor],
                                 ContextUtil.getConfiguration(job)).map(p => p._2).cache()
 
    // get phasing sets
    val keyedSides = sides.keyBy(_.getId).cache()
    val startAnch = an.keyBy(a => a.getEdge.getLeft)
      .join(keyedSides)
      .map(kv => kv._2._2.getPosition)
    val endAnch = an.keyBy(a => a.getEdge.getRight)
      .join(keyedSides)
      .map(kv => kv._2._2.getPosition)
    val phaseSets = startAnch.zip(endAnch)
      .map(r => (r, 1))
      .reduceByKey(_ + _)
      .groupBy(_._1._1.getContig)
      .flatMap(kv => {
        val (contig, ranges) = kv
        val builder = new PhaseSetBuilder(contig)
        
        ranges.map(t => {
          val ((s, e), c) = t

          (c, s.getBase, e.getBase)
        }).map(kv => (kv._1, kv._2.toLong, kv._3.toLong))
        .foreach(builder.add(_))

        builder.toSets()
      })
      .collect()

    // set contig and sample info in header
    SequenceGraphVCFOutputFormat.addSample(sample)
    
    sides.map(_.getPosition)
      .map(_.getContig)
      .map(_.toString)
      .distinct()
      .map(new java.lang.String(_))
      .collect()
      .foreach(SequenceGraphVCFOutputFormat.addContig(_))
    
    // get allele info
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

  def convertToVariantContext(position: Position,
                              alleles: List[AlleleGroup],
                              phaseSets: Array[PhaseSet],
                              sample: String): VariantContextWritable = {
    // convert alleles to Picard allele
    val alleleNoGroup = alleles.flatMap((r: AlleleGroup) => r.getAllele() match {
      case (a: SGAllele) => Some(a)
      case _ => None
    })

    val refBases = alleleNoGroup.head.getReferenceBases

    val alleleBases = alleleNoGroup.map(_.getBases())

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
    
    // check for phasing
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

    // make genotype context...
    val genotypeList: java.util.List[Genotype] = asList(List(genotype))
    val genotypeArrayList: java.util.ArrayList[Genotype] = new java.util.ArrayList(genotypeList)
    val sampleMap: java.util.Map[java.lang.String, java.lang.Integer] = asMap(Map(new java.lang.String(sample) -> new java.lang.Integer(0)))
    val sampleOrder: java.util.List[String] = asList(List(new java.lang.String(sample)))
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





















