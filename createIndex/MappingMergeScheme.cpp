#include <stdexcept>

#include <Log.hpp>
#include <util.hpp>

#include "MappingMergeScheme.hpp"

const size_t MappingMergeScheme::MAX_THREADS = 32;

MappingMergeScheme::MappingMergeScheme(const FMDIndex& index,
    const MappingScheme* mappingScheme, size_t genome): MergeScheme(index),
    index(index), genome(genome), mappingScheme(mappingScheme), threads(),
    queue(NULL), contigsToMerge(NULL) {
    
    // Nothing to do
    
}

MappingMergeScheme::~MappingMergeScheme() {
    
    join(); // Join our threads

    if(queue != NULL) {
        // Get rid of the queue if we made one.
        delete queue;
        queue = NULL;
    }
    
    if(contigsToMerge != NULL) {
        // And the collection of contigs to merge.
        delete contigsToMerge;
        contigsToMerge = NULL;
    }
    
}

ConcurrentQueue<Merge>& MappingMergeScheme::run() {
    
    if(queue != NULL) {
        // Don't let people call this twice.
        throw std::runtime_error(
            "Called run() twice on a MappingMergeScheme.");
    }
    
    // Grab the limits of the contig range belonging to this genome.
    auto genomeContigs = index.getGenomeContigs(genome);
    
    // Don't start more threads than we have contigs.
    size_t numThreads = std::min(MAX_THREADS, genomeContigs.second - 
        genomeContigs.first);
    
    Log::info() << "Running Mapping merge on " << numThreads << " threads" <<
        std::endl;
    
    // Make the queue of merges    
    queue = new ConcurrentQueue<Merge>(numThreads);
    
    // And one for the contigs to merge
    contigsToMerge = new ConcurrentQueue<size_t>(1);
    
    // Say we need to do every contig in this genome.
    
    for(size_t contig = genomeContigs.first; contig < genomeContigs.second; 
        contig++) {
      
        // Put each contig in the genome into the queue of work to do.
        // TODO: Can we just use an imaginary queue of numbers in a range?
        auto lock = contigsToMerge->lock();
        contigsToMerge->enqueue(contig, lock);
    }
    
    // Say we are done writing to the queue. Good thing it doesn't have a max
    // count. Otherwise we'd need a thread to dump in all the contigs.
    auto lock = contigsToMerge->lock();
    contigsToMerge->close(lock);

    // Start up a reasonable number of threads to do the work.
    for(size_t threadID = 0; threadID < numThreads; threadID++) {
        // All the mapTypes can use the same method here, since there's lots of
        // work scheduling and setup logic in common.
        
        // Start up a thread.
        threads.push_back(Thread(&MappingMergeScheme::generateMerges,
            this, contigsToMerge));
    }

    // Return a reference to the queue of merges, for our caller to do something
    // with.
    return *queue;
    
}

void MappingMergeScheme::join() {

    for(auto i = threads.begin(); 
        i != threads.end(); ++i) {
    
        // Join each thread.
        (*i).join();    
        
    }
    
    // Probably not a good idea to join the same threads twice.
    threads.clear();
    
    if(contigsToMerge != NULL) {
        auto lock = contigsToMerge->lock();
        if(!contigsToMerge->isEmpty(lock)) {
            // Don't hold a lock and throw an exception.
            lock.unlock();
        
            // Complain our thread pool somehow broke.
            throw std::runtime_error(
                "Jobs left in queue after threads have terminated.");
        } else {
            // Everything is fine, but we still need to unlock.
            lock.unlock();
        }
    }
    
}

void MappingMergeScheme::generateMerge(size_t queryContig, size_t queryBase, 
    size_t referenceContig, size_t referenceBase, bool orientation) const {
        
    // Where in the query do we want to come from? Always on the forward strand.
    // Correct offset to 0-based.
    TextPosition queryPos(queryContig * 2, queryBase - 1);
    
    // Where in the reference do we want to go to? May b e on the forward or
    // reverse strand. Correct text and offset for orientation, and correct
    // offset to 0-based (at the same time).
    TextPosition referencePos(referenceContig * 2 + orientation, 
        orientation ? (index.getContigLength(referenceContig) - referenceBase) :
        (referenceBase - 1));
    
    // Make a Merge between these positions.
    
    if(queryBase > index.getContigLength(queryContig) || queryBase == 0 ||
        referenceBase > index.getContigLength(referenceContig) || 
        referenceBase == 0) {
        
        // We're trying to merge something out of bounds.
        
        Log::critical() << 
            "****Tried to merge (to) an off-contig position!****" << 
            std::endl;
            
        throw std::runtime_error("Tried to merge (to) an off-contig position!");
    
    } else {
    
        Merge merge(queryPos, referencePos);
            
        // Send that merge to the queue.
        // Lock the queue.
        auto lock = queue->lock();
        // Spend our lock to add something to it.
        queue->enqueue(merge, lock);
        
    }
    
}

void MappingMergeScheme::generateMerges(
    ConcurrentQueue<size_t>* contigs) const {
    
    // Wait for a contig, or for there to be no more.
    auto contigLock = contigs->waitForNonemptyOrEnd();
    
    while(!contigs->isEmpty(contigLock)) {
        // We got a contig to do. Dequeue it and unlock.
        size_t queryContig = contigs->dequeue(contigLock);
        
        generateSomeMerges(queryContig);
        
        // Now wait for a new task, or for there to be no more contigs.
        contigLock = contigs->waitForNonemptyOrEnd();
        
    }
        
    // If we get here, there were no more contigs in the queue, and no
    // writers writing. Unlock it and finish the thread. TODO: Should the
    // queue just unlock if it happens to be empty?
    contigLock.unlock();
    
    // Close the output queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
}

void MappingMergeScheme::generateSomeMerges(size_t queryContig) const {
    
    
    // Set our task name to something descriptive.
    std::string taskName = "T" + std::to_string(genome) + "." + 
        std::to_string(queryContig);
        
    // Grab the contig as a string
    std::string contig = index.displayContig(queryContig);
    
    // How many bases are we trying to mapping?
    Log::info() << taskName << " mapping " << contig.size() << " bases." <<
        std::endl;
    
    // How many bases will we map?
    size_t mappedBases = 0;    
    
    // Go make all the mappings
    mappingScheme->map(contig, [&](size_t base, TextPosition mappedTo) {
        // Count each mapping
        mappedBases++;
        
        // Make the actual merge. Remember that position arguments need to be
        // 1-based.
        generateMerge(queryContig, base + 1, mappedTo.getContigNumber(),
            index.getContigOffset(mappedTo), mappedTo.getStrand());
    });
    
    Log::info() << taskName << " mapped " << mappedBases << "/" << 
        contig.size() << " bases." << std::endl;
}







