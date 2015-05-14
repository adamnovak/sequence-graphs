#ifndef CREDITSTRATEGY_HPP
#define CREDITSTRATEGY_HPP

#include "FMDPosition.hpp"
#include "FMDIndexView.hpp"
#include "FMDPositionGroup.hpp"

/**
 * A tool for applying credit between mapped bases, against a graph reference.
 */
class CreditStrategy {

public:
    
    /**
     * Make a new CreditStrategy using the given view of an FMDIndex. Can
     * optionally start disabled.
     */
    CreditStrategy(const FMDIndexView& view, bool enabled = true);
    
    /**
     * Apply credit for mapping the given string, by updating mappings in the
     * given vector of mappings.
     */
    void applyCredit(const std::string& query, 
        std::vector<Mapping>& toUpdate) const;
    
    /**
     * Apply credit for mapping the given string, by updating mappings in the
     * given vector of mappings between the specified mapped positions flanking
     * an unmapped region.
     */
    void applyCreditBetween(const std::string& query,
        std::vector<Mapping>& toUpdate, size_t leftAnchor,
        size_t rightAnchor) const;
        
    /**
     * How many mismatches will we tolerate when applying credit?
     */
    size_t maxMismatches = 2;
    
    /**
     * Set to false to not actually do anything.
     */
    bool enabled = true;

protected:

    /**
     * Run a right to left breadth-first greedy search for credit. Start at the
     * given index in the query, which is mapped to the given TextPosition in
     * the graph, and search to the left until maxDepth characters are searched,
     * or you run out of results. For every unique mapping between a query index
     * and a TextPosition, call the callback.
     */
    void breadthFirstSearch(const std::string& query, size_t queryStart,
        const TextPosition& referenceStart, size_t maxDepth, 
        std::function<void(size_t, TextPosition)> callback) const;

    /**
     * Keep the view of the index that we are working with.
     */
    const FMDIndexView& view;

};

#endif
