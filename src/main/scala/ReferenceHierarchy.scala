package edu.ucsc.genome

/**
 * Represents a reference hierarchy composed of reference structures at
 * different levels. Defined by a bottom-level index, a merging scheme for
 * building each level, and a possibly finished, or possibly null, graph of
 * Sides, Sites, Adjacencies, and Generalizations (with no negative vertex IDs).
 */
class ReferenceHierarchy(index: FMDIndex) {
    
    // We keep an array of levels, starting with our bottom-most string level.
    var levels: Array[ReferenceStructure] = 
        Seq(new StringReferenceStructure(index)).toArray
    
    /**
     * Load a ReferenceHierarchy from the given path, using the given
     * SparkContext. The saved hierarchy must have been created by the
     * createIndex C++ program.
     */
    def this(path: String) = {
        // Load an FMDIndex for a basename in that directory, and create the
        // bottom level.
        this(new FMDIndex(path + "/index.basename"))
        
        // Load the merged level 1, on the same index.
        levels :+= new MergedReferenceStructure(index, path + "/level1")
        
    }
    
    /**
     * Map the given string on the given level. Discards ambiguous mappings.
     */
    def map(level: Int, context: String): Seq[Option[Side]] = {
        levels(level).map(context)
    }
    
    /**
     * Map the given string on the given level, on the given side of each base.
     */
    def map(level: Int, context: String, face: Face): Seq[Option[Side]] = {
        levels(level).map(context, face)
    }
    
    /**
     * Map the given string to all levels. Discards ambiguous mappings. Returns
     * a sequence of sequences of mappings (Side or None).
     */
    def map(context: String): Seq[Seq[Option[Side]]] = {
        // Run the normal map for every level.
        (0 until levels.size).map(this.map(_, context))
    }
    
    /**
     * Save this ReferenceHierarchy to the specified directory.
     */
    def save(path: String) = {
        
        
    }
    
    
    
    
    
}


