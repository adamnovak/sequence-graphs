package edu.ucsc.genome

object SequenceGraphKryoProperties {

  def setupContextProperties() = {
    System.setProperty("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
    System.setProperty("spark.kryo.registrator", "edu.ucsc.genome.SequenceGraphKryoRegistrator")
    System.setProperty("spark.kryoserializer.buffer.mb", "4")
    System.setProperty("spark.kryo.referenceTracking", "false")
  }
}
