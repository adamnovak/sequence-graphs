#include "StatTracker.hpp"

#include <fstream>

StatTracker& StatTracker::operator+=(const StatTracker& other){
    std::lock(statsMutex, other.statsMutex);
    for(auto& pair : other.stats) {
        // For each stat in the other StatTracker, add it into our stats. We
        // don't just forward along to add because we would need a recursive
        // mutex.
        if(stats.count(pair.first)) {
            // We already have something for this stat, so add to it.
            stats[pair.first] += pair.second;
        } else {
            // We had nothing, so start out with this value.
            stats[pair.first] = pair.second;
        }
    }
    statsMutex.unlock();
    other.statsMutex.unlock();
}

void StatTracker::add(const std::string& stat, size_t amount) {
    // Lock stats so nobody else can use it.
    statsMutex.lock();
    if(stats.count(stat)) {
        // We already have something for this stat, so add to it.
        stats[stat] += amount;
    } else {
        // We had nothing, so start out with this value.
        stats[stat] = amount;
    }
    statsMutex.unlock();
}

StatTracker::StatTracker(const StatTracker& other) {
    // Lock both objects
    std::lock(statsMutex, other.statsMutex);
    // Copy over the stats
    stats = other.stats;  
    // Unlock  
    statsMutex.unlock();
    other.statsMutex.unlock();
}

StatTracker& StatTracker::operator=(const StatTracker& other) {
    // Lock both objects
    std::lock(statsMutex, other.statsMutex);
    // Copy over the stats
    stats = other.stats;  
    // Unlock  
    statsMutex.unlock();
    other.statsMutex.unlock();
    // Return reference
    return *this;
}

const size_t StatTracker::operator[](const std::string stat) const {
    statsMutex.lock();
    size_t toReturn = 0;
    if(stats.count(stat)) {
        // We already have something for this stat, so return that. We need to
        // use at since the map is const.
        toReturn = stats.at(stat);
    }

    // Unlock and send back the answer.
    statsMutex.unlock();
    return toReturn;
}

void StatTracker::save(const std::string& filename) const {
    // Open the output file
    std::ofstream file(filename.c_str());

    // Lock the stats so we get a consistent view.
    statsMutex.lock();
    for(auto& pair : stats) {
        // For each stat, save the name and value
        file << pair.first << "\t" << pair.second;
    }
    // Unlock the stats in case someone wants them later.
    statsMutex.unlock();

    // Close the file up and sync it to disk
    file.close();
}
