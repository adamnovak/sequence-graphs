package edu.ucsc.genome

import com.esotericsoftware.kryo.Kryo
import org.apache.spark.serializer.KryoRegistrator

class SequenceGraphKryoRegistrator extends KryoRegistrator {
  override def registerClasses(kryo: Kryo) {
    kryo.register(classOf[Adjacency], new AvroSerializer[Adjacency]())
    kryo.register(classOf[Anchor], new AvroSerializer[Anchor]())
    kryo.register(classOf[AlleleGroup], new AvroSerializer[AlleleGroup]())
    kryo.register(classOf[Side], new AvroSerializer[Side]())
  }
}
