// package.scala: contains implicit definitions that need to be in scope
// throuought the entire package.
package edu.ucsc

/**
 * Package object that contains implicits needed throughout the entire package,
 * and/or needed to use the package properly.
 *
 * Mostly this is orderings for Avro records. Avro defines a way to order
 * records, and indeed Avro-generated Java code all inherits from a base class,
 * comparable with itself, that actually can sort and compare things. However,
 * in Scala, this results in all Avro objects being comparable against the base
 * record type but not its subtypes; this results in the implicit production of
 * an Ordering for the base type but no Orderings for the subtypes.
 *
 * TODO: We should move to Scala 2.10+ and use macros to generate all of the
 * Avro orderings.
 *
 */
package object genome {
    
    /**
     * Define an ordering for Positions implicitly, rather than trying to edit the code-generated Avro code to add it.
     */
    implicit val positionOrdering = Ordering.by { (position: Position) =>
        (position.contig, position.base, position.face)
    }
    
    /**
     * Define an ordering for Sides: sort them by position.
     * TODO: Move this to Avro.
     */ 
    implicit val sideOrdering = Ordering.by { (side: Side) =>
        side.position
    }
    
}

