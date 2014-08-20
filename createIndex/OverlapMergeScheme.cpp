#include "OverlapMergeScheme.hpp"

#include <Log.hpp>

OverlapMergeScheme::OverlapMergeScheme(const FMDIndex& index,
    size_t minContext): MergeScheme(index), threads(), queue(NULL), 
    minContext(minContext) {
    
    // Nothing to do
    
}

OverlapMergeScheme::~OverlapMergeScheme() {
    
    join(); // Join our threads

    if(queue != NULL) {
        // Get rid of the queue if we made one.
        delete queue;
        queue = NULL;
    }
    
}

ConcurrentQueue<Merge>& OverlapMergeScheme::run() {
    
    if(queue != NULL) {
        // Don't let people call this twice.
        throw std::runtime_error(
            "Called run() twice on a OverlapMergeScheme.");
    }
    
    // Figure out how many threads to run. Let's do one per pair of genomes (not
    // counting pairs of the same genome), because more threads is better.
    
    // TODO: If we scale to large numbers of genomes, this may make OS thread
    // limits angry, and we may have to sequence these somehow so we don't try
    // to run them all at once.
    size_t threadCount = index.getNumberOfGenomes() * 
        (index.getNumberOfGenomes() - 1);
    
    Log::info() << "Running Overlap merge on " << threadCount << " threads" <<
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
            threads.push_back(Thread(&OverlapMergeScheme::generateMerges,
                this, i, j));
        }
    }

    // Return a reference to it.
    return *queue;
    
}

void OverlapMergeScheme::join() {

    for(auto i = threads.begin(); 
        i != threads.end(); ++i) {
    
        // Join each thread.
        (*i).join();    
        
    }
    
    // Probably not a good idea to join the same threads twice.
    threads.clear();
    
}

void OverlapMergeScheme::generateMerges(size_t targetGenome,
    size_t queryGenome) {
    
    // What's our thread name?
    std::string threadName = "T" + std::to_string(queryGenome) + "->" + 
        std::to_string(targetGenome);
    
    // Now we just have to generate some merges by mapping each contig in the
    // query to the target, and write them to the queue.
    
    // Keep track of total bases mapped and unmapped
    size_t basesMapped = 0;
    size_t basesUnmapped = 0;
    
    // Get the range of contigs we need to map.
    std::pair<size_t, size_t> contigRange = index.getGenomeContigs(queryGenome);
    
    // How many are in there?
    size_t contigsToDo = contigRange.second - contigRange.first;
    
    // Say how many we are doing.
    Log::info() << threadName << " has to map " << contigsToDo << " contigs" <<
        std::endl;
    
    for(size_t i = contigRange.first; i < contigRange.second; i++) {
        Log::debug() << threadName << " mapping contig " << i << std::endl;
        
        // Keep track of mapped and unmapped bases per contig.
        size_t basesMappedInContig = 0;
        size_t basesUnmappedInContig = 0;
        
        // Grab each contig as a string
        std::string contig = index.displayContig(i);
        
        // Map it to the target genome in both orientations, with a minimum
        // context as specified when we were constructed, and disambiguate.
        std::vector<Mapping> mappings = index.mapBoth(contig, targetGenome,
            minContext);
        
        for(size_t base = 0; base < mappings.size(); base++) {
            // For each base that we tried to map
            
            if(!mappings[base].is_mapped) {
                // Skip the unmapped ones
                basesUnmappedInContig++;
                continue;
            }
            
            // If we get here, we mapped a base.
            basesMappedInContig++;
            Log::debug() << threadName << " mapped base " << base << std::endl;
            
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
        
        // Put bases from this contig into total stats
        basesMapped += basesMappedInContig;
        basesUnmapped += basesUnmappedInContig;
        
        // How much is left to do?
        size_t contigsDone = i + 1 - contigRange.first;
        
        
        Log::info() << threadName << " mapped contig " << i << " (" << 
            basesMappedInContig << "|" << basesUnmappedInContig << ") " << 
            contigsDone << "/" << contigsToDo << std::endl;
        
    }
    
    // Close the queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
    // Report that we're done and how many bases we mapped.
    Log::info() << threadName << " finished (" << basesMapped << "|" << 
        basesUnmapped << ")" << std::endl;
}










