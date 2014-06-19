#include "SymmetricMergeScheme.hpp"

#include <Log.hpp>

SymmetricMergeScheme::SymmetricMergeScheme(const FMDIndex& index): 
    MergeScheme(index), threads(), queue(NULL) {
    
    // Nothing to do
    
}

SymmetricMergeScheme::~SymmetricMergeScheme() {
    
    join(); // Join our threads

    if(queue != NULL) {
        // Get rid of the queue if we made one.
        delete queue;
        queue = NULL;
    }
    
}

ConcurrentQueue<Merge>& SymmetricMergeScheme::run() {
    
    if(queue != NULL) {
        // Don't let people call this twice.
        throw std::runtime_error(
            "Called run() twice on a SymmetricMergeScheme.");
    }
    
    // Figure out how many threads to run. Let's do one per pair of genomes (not
    // counting pairs of the same genome), because more threads is better.
    
    // TODO: If we scale to large numbers of genomes, this may make OS thread
    // limits angry, and we may have to sequence these somehow so we don't try
    // to run them all at once.
    size_t threadCount = index.getNumberOfGenomes() * 
        (index.getNumberOfGenomes() - 1);
    
    Log::info() << "Running symmetric merge on " << threadCount << " threads" <<
        std::endl;
    
    // Make the queue    
    queue = new ConcurrentQueue<Merge>(threadCount);
    
    // Start some threads.
    for(size_t i = 0; i < index.getNumberOfGenomes(); i++) {
        for(size_t j = 0; j < index.getNumberOfGenomes(); j++) {
            if(i == j) {
                // Don't map genomes to themselves; that's silly.
                continue;
            }
            
            // Otherwise, for each pair of genomes in each order, start a thread
            // to map the one to the other.
            threads.push_back(std::thread(&SymmetricMergeScheme::generateMerges,
                this, i, j));
        }
    }

    // Return a reference to it.
    return *queue;
    
}

void SymmetricMergeScheme::join() {

    for(std::vector<std::thread>::iterator i = threads.begin(); 
        i != threads.end(); ++i) {
    
        // Join each thread.
        (*i).join();    
        
    }
    
    // Probably not a good idea to join the same threads twice.
    threads.clear();
    
}

void SymmetricMergeScheme::generateMerges(size_t targetGenome,
    size_t queryGenome) {
    
    // Now we just have to generate some merges by mapping each contig in the
    // query to the target, and write them to the queue.
    
    // Get the range of contigs we need to map.
    std::pair<size_t, size_t> contigRange = index.getGenomeContigs(queryGenome);
    
    for(size_t i = contigRange.first; i < contigRange.second; i++) {
        Log::debug() << "Merge thread mapping contig " << i << " to genome " << 
            targetGenome << std::endl;
        
        // Grab each contig as a string
        std::string contig = index.displayContig(i);
        
        // Map it to the target genome in both orientations, and disambiguate.
        std::vector<Mapping> mappings = index.mapBoth(contig, targetGenome);
        
        
        for(size_t base = 0; base < mappings.size(); base++) {
            // For each base that we tried to map
            
            if(!mappings[base].is_mapped) {
                // Skip the unmapped ones
                continue;
            }
            
            Log::debug() << "Mapped base " << base << std::endl;
            
            // Produce a merge between the base we're looking at on the forward
            // strand of this contig, and the location (and strand) it mapped to
            // in the other genome.
            Merge merge(TextPosition(i * 2, base), mappings[base].location);
            
            // Send that merge to the queue.
            // Lock the queue.
            auto lock = queue->lock();
            // Spend our lock to add something to it.
            queue->enqueue(merge, lock);
        }
        
    }
    
    // Close the queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
}










