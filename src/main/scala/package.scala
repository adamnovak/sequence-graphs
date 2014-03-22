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
     * Define an ordering for Positions implicitly, rather than trying to edit
     * the code-generated Avro code to add it.
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
    
    /**
     * Allow a Side to be used in place of its Long ID. We do this as a down-
     * conversion (rather than converting Longs up to Sides) because otherwise
     * we'd have to make up info for the other Side fields (or at least null
     * them out), and could get into trouble. API users will just have to figure
     * out that they can use Sides wherever a Long is called for.
     */
    implicit def Side2Long(side: Side) = side.id
    
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
     * Magically wrap Sites in HasEdge objects.
     *
     * TODO: Upgrade to Scala 2.10 Implicit Classes(?) instead of this. But we'd
     * have to get rid of case class-ness.
     */
    implicit def Site2HasEdge(site: Site): HasEdge = {
        SiteEdge(site)
    }
    
    /**
     * Magically wrap Breakpoints in HasEdge objects.
     *
     * TODO: Upgrade to Scala 2.10 Implicit Classes(?) instead of this. But we'd
     * have to get rid of case class-ness.
     */
    implicit def Breakpoint2HasEdge(breakpoint: Breakpoint): HasEdge = {
        BreakpointEdge(breakpoint)
    }
    
    /**
     * Magically wrap Generalization in HasEdge objects.
     *
     * TODO: Upgrade to Scala 2.10 Implicit Classes(?) instead of this. But we'd
     * have to get rid of case class-ness.
     */
    implicit def Generalization2HasEdge(
        generalization: Generalization): HasEdge = {
        
        GeneralizationEdge(generalization)
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
            string.reverse.map {
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
     * Magically allow us to do nice things with Positions.
     */
    implicit class RichPosition(position: Position) {
        /** 
         * Flip Position faces with unary !
         */
        def unary_! : Position = {
            // Flip the face around
            val newFace = position.face match {
                case Face.LEFT => Face.RIGHT
                case Face.RIGHT => Face.LEFT
            }
            
            // Make a new Position
            new Position(position.contig, position.base, newFace)
        }
    }
}

