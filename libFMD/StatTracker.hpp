#ifndef STATTRACKER_HPP
#define STATTRACKER_HPP
// StatTracker.hpp: Defines a class for multi-threaded mapping scheme stat
// tracking

#include <map>
#include <string>
#include <mutex>

/**
 * A MappingScheme can have an associated StatTracker, for tracking things like
 * the number of aligned or unaligned bases. A StatTracker tracks one or more
 * integer stats (defined by string names), with atomic update methods, and can
 * generate a report at the end of the program's execution which can then be
 * merged together across runs.
 */
class StatTracker {
public:
    // We need our own copy constructors and assignment operators because we
    // have a mutex.
    
    /**
     * Default default constructor should exist.
     */
    StatTracker() = default;
    
    /**
     * Copy constructor needs to lock the source before copying.
     */
    StatTracker(const StatTracker& other);
    
    /**
     * Assignment operator needs to lock both source and destination.
     */
    StatTracker& operator=(const StatTracker& other);
    
    /**
     * Add another StatTracker into this one in a thread-safe but not
     * necessarily atomic way.
     */
    StatTracker& operator+=(const StatTracker& other);
    
    /**
     * Atomically add an amount to a stat.
     */
    void add(const std::string& stat, size_t amount);
    
    /**
     * Retrieve the value of a stat.
     */
    const size_t operator[](const std::string stat) const;
    
    /**
     * Save a TSV-format report of <stat>\t<value> pairs to the given file.
     */
    void save(const std::string& filename) const;
        
private:
    /**
     * Holds all the stat values by name.
     */
    std::map<std::string, size_t> stats;
    
    /**
     * Control access to the stats container. Note that this is non-recursive!
     * Needs to be mutable so we can lock a const thing to read from it safely.
     */
    mutable std::mutex statsMutex;
};

#endif
