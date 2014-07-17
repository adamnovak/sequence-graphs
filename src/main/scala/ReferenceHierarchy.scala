package edu.ucsc.genome
import org.ga4gh.FMDIndex

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
     * Map the given string on the given level. Only keeps mappings with the
     * given minimum context or longer. Discards ambiguous mappings.
     */
    def mapLevel(level: Int, context: String, minContext: Int = 0): 
        Seq[Option[Side]] = {
        
        levels(level).map(context, minContext)
    }
    
    /**
     * Map the given string on the given level, on the given face of each base.
     */
    def mapFaceLevel(level: Int, context: String, face: Face,
        minContext: Int = 0): Seq[Option[Side]] = {

        levels(level).mapFace(context, face, minContext)
    }
    
    /**
     * Map the given string to all levels. Discards ambiguous mappings. Returns
     * a sequence of sequences of mappings (Side or None). Only keeps mappings
     * with the given minimum context or longer.
     */
    def map(context: String, minContext: Int = 0): Seq[Seq[Option[Side]]] = {
        // Run the normal map for every level.
        (0 until levels.size).map(this.mapLevel(_, context, minContext))
    }
    
    /**
     * Map the given string to all levels, on the given face of each base.
     * Returns a sequence of sequences of mappings (Side or None). Only keeps
     * mappings with the given minimum context or longer.
     */
    def mapFace(context: String, face: Face, minContext: Int = 0):
        Seq[Seq[Option[Side]]] = {
        
        // Run the normal map for every level.
        (0 until levels.size).map(this.mapFaceLevel(_, context, face, 
            minContext))
    }
}


