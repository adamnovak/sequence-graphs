// package.scala: contains implicit definitions that need to be in scope
// throuought the entire package.
package edu.ucsc

/**
 * Package object that contains implicits needed throughout the entire package,
 * and/or needed to use the package properly.
 *
 * Mostly this is:
 *
 * Orderings for Avro records. Avro defines a way to order
 * records, and indeed Avro-generated Java code all inherits from a base class,
 * comparable with itself, that actually can sort and compare things. However,
 * in Scala, this results in all Avro objects being comparable against the base
 * record type but not its subtypes; this results in the implicit production of
 * an Ordering for the base type but no Orderings for the subtypes.
 *
 * Implicit constructors for Avro records. We don't get to add them ourselves to
 * the classes, so we make do here and provide some syntactic sugar.
 *
 * TODO: We should move to Scala 2.10+ and use macros to generate all of the
 * Avro orderings.
 *
 */
package object genome {
    
    /**
     * Define an ordering for Sides: sort them by position.
     * TODO: Move this to Avro.
     */ 
    implicit val sideOrdering = Ordering.by { (side: Side) =>
        (side.coordinate, side.face)
    }
    
    /**
     * Auto-Magically make PloidyBounds for fixed ploidies.
     */
    implicit def Int2PloidyBounds(ploidy: Int) = {
        new PloidyBounds(ploidy, null, null)
    }
    
    /**
     * Auto-Magically make PloidyBounds for 2-tuple ranges.
     */
    implicit def IntTuple2PloidyBounds(bounds: (Int, Int)) = {
        new PloidyBounds(bounds._1, bounds._2, null)
    }
    
    /**
     * Magically wrap AlleleGroups in HasEdge objects.
     *
     * TODO: Upgrade to Scala 2.10 Implicit Classes(?) instead of this. But we'd
     * have to get rid of case class-ness.
     */
    implicit def AlleleGroup2HasEdge(alleleGroup: AlleleGroup): HasEdge = {
        AlleleGroupEdge(alleleGroup)
    }
    
    /**
     * Magically wrap Adjacencies in HasEdge objects.
     *
     * TODO: Upgrade to Scala 2.10 Implicit Classes(?) instead of this. But we'd
     * have to get rid of case class-ness.
     */
    implicit def Adjacency2HasEdge(adjacency: Adjacency): HasEdge = {
        AdjacencyEdge(adjacency)
    }
    
    /**
     * Magically wrap Anchors in HasEdge objects.
     *
     * TODO: Upgrade to Scala 2.10 Implicit Classes(?) instead of this. But we'd
     * have to get rid of case class-ness.
     */
    implicit def Anchor2HasEdge(anchor: Anchor): HasEdge = {
        AnchorEdge(anchor)
    }
    
    /**
     * Magically wrap Blocks in HasEdge objects.
     *
     * TODO: Upgrade to Scala 2.10 Implicit Classes(?) instead of this. But we'd
     * have to get rid of case class-ness.
     */
    implicit def Block2HasEdge(block: Block): HasEdge = {
        BlockEdge(block)
    }
    
    /**
     * Magically provide DNA methods on strings.
     */
    implicit class DNAString(string: String) {
        /**
         * Return the reverse complement of this promoted string. All letters in
         * the string must be {A, C, G, T, N}, with no gaps or lower-case
         * characters.
         */
        def reverseComplement: String = {
            string.reverse.map(_.reverseComplement)
        }
    }
    
    /**
     * Magically provide DNA methods on characters.
     */
    implicit class DNAChar(char: Char) {
        /**
         * Return the reverse complement of this promoted char. The character
         * must be {A, C, G, T, N}, with no gaps or lower-case characters.
         */
        def reverseComplement: Char = {
            char match {
                case 'A' => 'T'
                case 'C' => 'G'
                case 'G' => 'C'
                case 'T' => 'A'
                case 'N' => 'N'
            }
        }
    }
    
    /**
     * Magically allow us to do nice things with Faces.
     */
    implicit class RichFace(face: Face) {
        /** 
         * Flip Faces with unary !
         */
        def unary_! : Face = face match {
            case Face.LEFT => Face.RIGHT
            case Face.RIGHT => Face.LEFT
        }
    }
    
    /**
     * Magically allow us to do nice things with Sides.
     */
    implicit class RichSide(side: Side) {
        /** 
         * Flip Side faces with unary !
         */
        def unary_! : Side = {
            // Flip the face around
            val newFace = side.face match {
                case Face.LEFT => Face.RIGHT
                case Face.RIGHT => Face.LEFT
            }
            
            // Make a new Side
            new Side(side.coordinate, newFace)
        }
    }
}

