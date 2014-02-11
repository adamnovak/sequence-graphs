package edu.ucsc.genome.ExportVCF

import fi.tkk.ics.hadoop.bam.{KeyIgnoringVCFOutputFormat, VCFFormat}
import org.broadinstitute.variant.vcf.{VCFHeader, VCFFormatHeaderLine, VCFInfoHeaderLine, VCFContigHeaderLine, VCFConstants, VCFHeaderLine, VCFHeaderLineType}
import scala.collection.JavaConversions._

object SequenceGraphVCFOutputFormat {
  
  var samples = Set[String]()
  var contigs = List[String]()
  
  def addSample (sample: String) {
    samples = samples + sample
  }
  
  def addContig (contig: String) {
    contigs = contig :: contigs
  }

}


class SequenceGraphVCFOutputFormat[K]
  extends KeyIgnoringVCFOutputFormat[K](VCFFormat.valueOf("VCF")) {
  
  // make header
  var headerLines = Set[VCFHeaderLine](new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"),
                                       new VCFInfoHeaderLine(VCFConstants.PHASE_SET_KEY, 1, VCFHeaderLineType.Integer, "Phase set")
                                     ) /*++ SequenceGraphVCFOutputFormat.contigs.map(n => Map(n -> n))
    .zip(0 until SequenceGraphVCFOutputFormat.contigs.length)
    .map(pair => new VCFContigHeaderLine(asMap(pair._1), pair._2).asInstanceOf[VCFHeaderLine])
    .toSet
    * TODO: fix the code above-seems to be an issue in the VCFContigHeaderLine
    * */

  val hdr = new VCFHeader(setAsJavaSet(headerLines), setAsJavaSet(SequenceGraphVCFOutputFormat.samples.map(new java.lang.String(_))))

  setHeader(hdr)
}
