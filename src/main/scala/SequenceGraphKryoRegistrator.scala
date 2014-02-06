package edu.ucsc.genome

import com.esotericsoftware.kryo.Kryo
import org.apache.spark.serializer.KryoRegistrator

/**
 * Class to register Sequence Graph types for Kryo
 * serialization/deserialization.
 */
class SequenceGraphKryoRegistrator extends KryoRegistrator {
  override def registerClasses(kryo: Kryo) {
    // Register all the Sequence Graph parts
    kryo.register(classOf[Adjacency], new AvroSerializer[Adjacency]())
    kryo.register(classOf[Anchor], new AvroSerializer[Anchor]())
    kryo.register(classOf[AlleleGroup], new AvroSerializer[AlleleGroup]())
    kryo.register(classOf[Side], new AvroSerializer[Side]())
  }
}
