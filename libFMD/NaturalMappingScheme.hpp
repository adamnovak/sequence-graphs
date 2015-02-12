#ifndef NATURALMAPPINGSCHEME_HPP
#define NATURALMAPPINGSCHEME_HPP

#include "MappingScheme.hpp"
#include "Mapping.hpp"
#include "IntervalIndex.hpp"
#include "Log.hpp"
#include "Matching.hpp"

#include <iomanip>

#include <unordered_map>

/**
 * Mapping scheme implementing Benedict's "natural" mapping scheme. If all
 * the unique-in-the-reference strings overlapping a query base agree about
 * where the abse should be, and you can't get from the base off the end of
 * the query with multiple results, map the base.
 *
 * This version counts edit distance when linking maximal matchings.
 */
class NaturalMappingScheme: public MappingScheme {

public:
    /**
     * Inherit the constructor.
     */
    using MappingScheme::MappingScheme;
    
    /**
     * Map the given query string according to the natural mapping scheme with
     * optional inexact credit. When a mapping is found, the callback function
     * will be called with the query base index, and the TextPosition to which
     * it maps in the forward direction.
     *
     */
    virtual void map(const std::string& query,
        std::function<void(size_t, TextPosition)> callback) const override;
        
    
    // Now come the scheme parameters and their default values.
    
    /**
     * What's the minimum factor by which the maximum unique string length must
     * exceed the minimum?
     */
    double multContext = 0.0;
    
    /**
     * What's the minimum max matching length in bases to use? Doesn't interact
     * with synteny runs exactly right.
     */
    size_t minContext = 0;
    
    /**
     * Should maximal unique matchings not in an admissible run still produce
     * conflict?
     */
    bool conflictBelowThreshold = false;
    
    /**
     * How many mismatches can be tolerated when applying credit?
     */
    size_t z_max = 0;
    
    /**
     * Should credit be used to create more mappings?
     */
    bool credit = false;
    
    /**
     * What is the minimum length that a maximal unique match has to be to even
     * contribute to conflicts? Maximal unique matches shorter than this will be
     * completely ignored, and won't be able to stop any bases from mapping on
     * other maximal unique matches.
     */
    size_t ignoreMatchesBelow = 0;
    
    /**
     * How many non-overlapping minimal unique matches does a maximal unique
     * match (or region of a synteny block) need in order to map?
     *
     * This is a lower bound on the maximal unique match's Hamming distance from
     * every other place in the reference except the place it mapped to.
     */
    size_t minHammingBound = 0;
    
    /**
     * How many mismatches are we allowed to have in gaps between maximal unique
     * matches that are working together to map?
     */
    size_t maxHammingDistance = 0;
    
    /**
     * How big of an alignment will we tolerate doing for edit distance
     * calculations?
     */
    size_t maxAlignmentSize = 500;
    
    /**
     * Should we enable potentially unstable (as in, mappings can change on
     * extension, not as in may crash the program) mapping hacks to increase
     * coverage when mapping reads?
     */
    bool unstable = false;
    
protected:
    /**
     * Map the given query string, producing a vector of Mappings. Does not
     * include credit yet.
     */
    std::vector<Mapping> naturalMap(const std::string& query) const;
    
    /**
     * Find all of the maximal unique matchings between query string characters
     * and the reference, in descending order by left endpoint. Properly handles
     * merged reference positions.
     */
    std::vector<Matching> findMaxMatchings(const std::string& query) const;
        
    /**
     * Find all of the minimal unique matchings between query string characters
     * and the reference, in descending order by left endpoint. Properly handles
     * merged reference positions.
     */
    std::vector<Matching> findMinMatchings(const std::string& query) const;
    
    /**
     * Produce a graph from each MUM to the MUMs it connects to, with the
     * mismatch gap cost of the connection. Includes self edges at cost 0. Input
     * matchings must not contain each other and must be in ascending order of
     * start position. Note that this is the reverse order of the
     * mindMaxMatchings method! TODO: Typedef this return type.
     */
    std::unordered_map<Matching, std::vector<std::pair<Matching, size_t>>>
        generateMaxMatchingGraph(std::vector<Matching> maxMatchings,
        const std::string& query) const;
        
    /**
     * Given a graph from each MUM to the MUMs it connects to, this the mismatch
     * gap cost of the connection, invert all the directed edges in the graph,
     * creating a new graph. Edges from a node may be in any order.
     */
    std::unordered_map<Matching, std::vector<std::pair<Matching, size_t>>>
        invertGraph(const std::unordered_map<Matching,
        std::vector<std::pair<Matching, size_t>>>& graph) const;
        
    /**
     * Look up the max matching that contains a given min matching.
     */
    const Matching& getMaxMatching(const IntervalIndex<Matching>& maxMatchings,
        const Matching& minMatching) const;
    
    /**
     * Given an index of max matchings, and a vector of min matchings in
     * ascending order (opposite of findMinMatchings), produce a map from max
     * matchings to the indexed min matchings they contain, in ascending order.
     */
    std::unordered_map<Matching, IntervalIndex<Matching>> assignMinMatchings(
        const IntervalIndex<Matching>& maxMatchings,
        const std::vector<Matching>& minMatchings) const;
        
    /**
     * For each min matching, and for each cost value <= maxHammingDistace,
     * calculate the number of nonoverlapping minimal unique matchings that can
     * be chained together, going in a certain direction (forward or reverse).
     *
     * Depends on the graph of max matching connectivity (with edges in any
     * order), the assignments of min matchings to max matchings (in ascending
     * order), and the list of max matchings (in ascending order). The graph
     * should contain self edges at cost 0.
     */
    std::unordered_map<Matching, std::vector<size_t>> getMinMatchingChains(
        const std::unordered_map<Matching, 
        std::vector<std::pair<Matching, size_t>>>& maxMatchingGraph,
        const std::unordered_map<Matching, IntervalIndex<Matching>>&
        minsForMax, const std::vector<Matching>& maxMatchings,
        bool isForward) const;
        
    /**
     * Given vectors of max and min matchings in ascending order as well as the
     * query, produce a vector of only max matchings that have min matchings in
     * sufficiently good (>= minHammingBound) synteny runs.
     */
    std::set<Matching> maxMatchingsInSyntenyRuns(const std::vector<Matching>&
        maxMatchings, const std::vector<Matching>& minMatchings,
        const std::string& query) const;
    
    /**
     * Compute the edit distance between the specified region of the query and
     * the specified region of the reference. If the distance would be greater
     * than the threshold, may return the threshold value instead. Also returns
     * the threshold value if the alignment that would need to be done is larger
     * than maxAlignmentSize.
     *
     * If leftJustify is false, do not charge for extra reference bases on the
     * left. If rightJustify is false, do not charge for extra reference bases
     * on the right.
     */
    size_t countEdits(const std::string& query, size_t queryStart,
        size_t queryLength, TextPosition referenceStart, size_t referenceLength,
        int64_t threshold = -1, bool leftJustify = true,
        bool rightJustify = true) const;
    
    
};

#endif
