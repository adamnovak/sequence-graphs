#include "adjacencyComponentUtil.hpp"
#include <Log.hpp>

// Pull in VFLib for graph matching, to deduplicate isomorphic components.
#include <argraph.h>
#include <argedit.h>
#include <vf2_sub_state.h>

#include <algorithm>
#include <unordered_set>
#include <fstream>

std::vector<std::vector<stPinchEnd>> 
getAdjacencyComponents(
    stPinchThreadSet* threadSet
) {

    Log::info() << "Making adjacency component list..." << std::endl;

    // Make an empty vector of components to populate.
    std::vector<std::vector<stPinchEnd>> toReturn;
    
    // Get all the adjacency components.
    stList* adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(
        threadSet);
        
    // Get an iterator over them
    stListIterator* componentIterator = stList_getIterator(adjacencyComponents);
    
    // Grab the first component
    stList* component = (stList*) stList_getNext(componentIterator);
    
    while(component != NULL) {
        
        // Make a vector to hold the ends in the component.
        std::vector<stPinchEnd> ends;
        
        // Get an iterator over its contents
        stListIterator* endIterator = stList_getIterator(component);
        
        // Get the first pinch end in the component
        stPinchEnd* pinchEnd = (stPinchEnd*) stList_getNext(endIterator);
        
        while(pinchEnd != NULL) {
            // For each pinch end in the component
            
            Log::trace() << LOG_LAZY("Observed end " << pinchEnd << 
                " of block " << stPinchEnd_getBlock(pinchEnd) << 
                " of degree " << 
                stPinchBlock_getDegree(stPinchEnd_getBlock(pinchEnd)) << 
                std::endl);
            
            // Put it in the vector that represents the component. Make sure to
            // copy it since the actual ends pointed to here get destroyed when
            // the list we're iterating over does.
            ends.push_back(*pinchEnd);
        
            // Look at the next end
            pinchEnd = (stPinchEnd*) stList_getNext(endIterator);
        }
        
        // Clean up the iterator
        stList_destructIterator(endIterator);
        
        // Put this component in the vector of all components.
        toReturn.push_back(std::move(ends));
        
        // Look at the next component
        component = (stList*) stList_getNext(componentIterator);
    }
    
    // Clean up the iterator
    stList_destructIterator(componentIterator);
    
    // And the entire list while we're at it
    stList_destruct(adjacencyComponents);
    
    // Give back the converted data structure
    return toReturn;

}

std::map<size_t, size_t>
getAdjacencyComponentSpectrum(
    const std::vector<std::vector<stPinchEnd>>& components
) {
    
    Log::info() << "Making adjacency component spectrum..." << std::endl;

    // Make the map we're going to return
    std::map<size_t, size_t> toReturn;

    for(auto component : components) {
        
        // Get its size
        size_t componentSize = component.size();
        
        if(!toReturn.count(componentSize)) {
            // This is the first component of this size we have found
            toReturn[componentSize] = 1;
        } else {
            // We found another one
            toReturn[componentSize]++;
        }
    }
    
    // Give back the map
    return toReturn;
    
}

std::vector<std::vector<stPinchEnd>>
filterComponentsBySize(
    const std::vector<std::vector<stPinchEnd>>& components,
    size_t size,
    std::function<bool(size_t, size_t)> comparison
) {
    
    Log::info() << "Selecting size " << size << " components" << std::endl;

    // Make a vector to return
    std::vector<std::vector<stPinchEnd>> toReturn;
    
    std::copy_if(components.begin(), components.end(), 
        std::back_inserter(toReturn), [&](std::vector<stPinchEnd> v) {
            // Grab only the vectors that compare correctly to size (according
            // to the passed predicate, which defaults to equality.
            return comparison(v.size(), size);
        });
        
    return toReturn;
}

std::vector<std::vector<stPinchSegment*>>
getAllPaths(
    stPinchEnd* start,
    stPinchEnd* end
) {
    
    Log::debug() << LOG_LAZY("Getting all paths from block " << 
        stPinchEnd_getBlock(start) << " orientation " << 
        stPinchEnd_getOrientation(start) << " to block " << 
        stPinchEnd_getBlock(end) << " orientation " << 
        stPinchEnd_getOrientation(end) << std::endl);

    // make the vector of paths we are going to return.
    std::vector<std::vector<stPinchSegment*>> toReturn;
    
    // Iterate over the threads in the first end's block
    stPinchBlockIt segmentsIterator = stPinchBlock_getSegmentIterator(
        stPinchEnd_getBlock(start));
        
    // Grab the first segment
    stPinchSegment* startSegment = stPinchBlockIt_getNext(&segmentsIterator);
    
    while(startSegment != NULL) {
    
        // Start a traversal from that segment
        stPinchSegment* segment = startSegment;
        
        // Make a vector to include all the segments except the bookending ones.
        std::vector<stPinchSegment*> path;
        
        // Keep track of the direction we are going. 0 is forwards.
        bool direction = stPinchEnd_getOrientation(start);
        
        Log::trace() << LOG_LAZY("Starting direction " << direction <<
            " from segment " << startSegment <<  " in block " << 
            stPinchEnd_getBlock(start) << std::endl);
        
        while(true) {
            // While we haven't made it to the block of the other end
            
            if(!direction) {
                // Go forwards (towards 3')
                segment = stPinchSegment_get3Prime(segment);
            } else {
                // Go backwards (towards 5')
                segment = stPinchSegment_get5Prime(segment);
            }
            
            
            if(segment == NULL) {
                // We ran off the end of the thread without getting to the thing
                // we were supposed to hit. Don't keep this path, and try the
                // next start segment.
                Log::error() << "Path escaped!" << std::endl;
                break;
            }
            
            // Get the orientation of this segment. 0 is forwards.
            bool segmentOrientation = stPinchSegment_getBlockOrientation(
                segment);
            
            // We know segments are connected 5' to 3' along a whole thread, so
            // we never need to change direction.
            
            Log::trace() << LOG_LAZY("Visiting segment " << segment << 
                " orientation " << segmentOrientation << " in block " << 
                stPinchSegment_getBlock(segment) << std::endl);
            
            if(stPinchSegment_getBlock(segment) == end->block && 
                direction != stPinchEnd_getOrientation(end)) {
                // We hit the other end. For orientation, if direction is true,
                // we were going backwards and have hit the 5' end. If direction
                // is false, we were going forwards and have hit the 3' end. End
                // orientation is true if it's the 3' end. So we need direction
                // and end orientation to not match if we are to detect the
                // correct end.
                
                // This might not matter for the size 2 adjacency component case
                // because any block we hit will be the right one. But if we
                // ever use this in more complex cases we will care.
                
                // Keep our path.
                toReturn.push_back(path);
                
                Log::debug() << LOG_LAZY("Path finished successfully with " << 
                    path.size() << " segments" << std::endl);
                
                // Stop this path.
                break;
            }
            
            // If we didn't hit the other end yet, keep this new segment on our
            // path.
            path.push_back(segment);
            
        }
    
        // Try the next segment to get the path from it.
        startSegment = stPinchBlockIt_getNext(&segmentsIterator);
    }
    
    // Return the vector of paths.
    return toReturn;

}

std::vector<int64_t>
getIndelLengths(
    const std::vector<std::vector<stPinchEnd>>& components
) {

    Log::info() << "Getting indel lengths from " << components.size() << 
        " components..." << std::endl;

    // Make the vector of lengths we will populate.
    std::vector<int64_t> toReturn;
    
    for(auto component : components) {
        // Look at every component
        
        if(component.size() != 2) {
            // Complain we got a bad sized component in here
            throw std::runtime_error(
                std::string("Component has incorrect size ") + 
                std::to_string(component.size()) + " for an indel.");
        }
        
        if(stPinchBlock_getDegree(stPinchEnd_getBlock(&component[0])) != 2 ||
            stPinchBlock_getDegree(stPinchEnd_getBlock(&component[1])) != 2) {
            // If there aren't exactly two segments on both ends, this isn't a
            // nice simple indel/substitution. Skip it.
            
            Log::error() << "Skipping component; Degrees are " << 
                stPinchBlock_getDegree(stPinchEnd_getBlock(&component[0])) << 
                " and " << 
                stPinchBlock_getDegree(stPinchEnd_getBlock(&component[1])) << 
                std::endl;
            continue;
            
            // TODO: account for ends of contigs lining up across from things
            // and having no paths across.
        }
        
        // Get two paths of 0 or more segments
        std::vector<std::vector<stPinchSegment*>> paths = 
            getAllPaths(&component[0], &component[1]);
            
        if(paths.size() != 2) {
            // This isn't just an indel. It might be an indel within a
            // duplication, or an indel for one out of a set of merged genomes,
            // but without just 2 paths we can't really define indel length.
            // TODO: Handle the case where all paths but one are one length, and
            // one is another length.
            Log::error() << "Skipping component; got " << paths.size() << 
                " paths instead of 2" << std::endl;
            continue;
        }
        
        // Keep the length of the segments on each path. There will always be 2
        // paths.
        std::vector<int64_t> pathLengths;
        
        for(auto path : paths) {
            // For each path, we want to track the length
            int64_t pathLength = 0;
            for(size_t i = 0; i < path.size(); i++) {
                // Add in each segment
                pathLength += stPinchSegment_getLength(path[i]);
                
                if(stPinchSegment_getLength(path[i]) > 100000000000000) {
                    // I saw some pretty wrong indel lengths.
                    throw std::runtime_error("Segment stupidly long");
                }
            }
            
            Log::debug() << "Path length: " << pathLength << std::endl;
            
            // Record the length of this path.
            pathLengths.push_back(pathLength);
        }
        
        if(pathLengths[0] < 0 || pathLengths[1] < 0) {
            // Maybe wrong indel lengths are from negative path lengths?
            throw std::runtime_error("Negative length path");
        }
        
        if(pathLengths[0] > 100000000000000 || 
            pathLengths[1] > 100000000000000) {
            // Maybe wrong indel lengths are from huge path lengths?
            throw std::runtime_error("Path stupidly long");
        }
        
        if(pathLengths[0] > pathLengths[1]) {
            // This is the way the indel goes
            toReturn.push_back(pathLengths[0] - pathLengths[1]);
        } else {
            // It goes the other way around.
            toReturn.push_back(pathLengths[1] - pathLengths[0]);
        }
    }
    
    return toReturn;

}

size_t
countTandemDuplications(
    const std::vector<std::vector<stPinchEnd>>& components
) {

    Log::info() << "Counting tandem duplications in " << components.size() << 
        " components..." << std::endl;

    // How many have we found so far?
    size_t tandemDuplications = 0;

    for(auto component : components) {
        // Is this component a tandem duplication?
        bool isTandem = false;
    
        for(size_t i = 0; i < component.size(); i++) {
            // Go through all the ends
            stPinchEnd end1 = component[i];
            
            Log::debug() << "End " << stPinchEnd_getOrientation(&end1) << 
                " of block " << stPinchEnd_getBlock(&end1) << std::endl;
            
            for(size_t j = 0; j < i; j++) {
                // And all the other ends
                stPinchEnd end2 = component[j];
                
                if(stPinchEnd_getBlock(&end1) == stPinchEnd_getBlock(&end2)) {
                    // We have two ends that share a block.
                    
                    Log::debug() << "We have two ends of block " << 
                        stPinchEnd_getBlock(&end1) << std::endl;
                    
                    // Get the ends attached to end 1
                    stSet* connectedEnds = 
                        stPinchEnd_getConnectedPinchEnds(&end1);
                        
                    // These ends are in there by address, so we have to scan
                    // for ours.
                    
                    // Get an iterator over the set.
                    stSetIterator* iterator = stSet_getIterator(connectedEnds);
                    
                    stPinchEnd* other = (stPinchEnd*) stSet_getNext(iterator);
                    while(other != NULL) {
                        // Go through all the things attached to end1
                        
                        Log::debug() << "Connection to block " << 
                            stPinchEnd_getBlock(other) << " end " << 
                            stPinchEnd_getOrientation(other) << std::endl;
                        
                        if(stPinchEnd_getBlock(other) == 
                            stPinchEnd_getBlock(&end2) && 
                            stPinchEnd_getOrientation(other) == 
                            stPinchEnd_getOrientation(&end2)) {
                            
                            Log::debug() << "...which counts!" << std::endl;
                            
                            // This end that the first end is connected to looks
                            // exactly like the second end. Call this a tandem
                            // duplication.
                            isTandem = true;
                            
                            // TODO: Break out of like 3 loops now, so we don't
                            // check all the other end pairs or somehow call two
                            // tandem duplications in one component.
                            
                        } else {
                            Log::debug() << "...which isn't block " << 
                                stPinchEnd_getBlock(&end2) << " end " << 
                                stPinchEnd_getOrientation(&end2) << std::endl;
                        }
                                                
                        other = (stPinchEnd*) stSet_getNext(iterator);
                    }
                    
                    // Clean up the iterator
                    stSet_destructIterator(iterator);
                        
                    // Clean up our connected ends set.
                    stSet_destruct(connectedEnds);
                }
            }
        }
        
        // Now we know if the component is a tandem duplication or not. Add to
        // the sum if applicable.
        tandemDuplications += isTandem;
        
        if(!isTandem) {
            Log::output() << "Non-tandem-duplication: " << std::endl;
            
            for(size_t i = 0; i < component.size(); i++) {
                // Go through all the ends
                stPinchEnd end1 = component[i];
                
                // Get the block for each
                stPinchBlock* block = stPinchEnd_getBlock(&end1);
                
                // And a segment for the block
                stPinchSegment* segment = stPinchBlock_getFirst(block);
                
                // Report the block (already done 1-based I think?)
                Log::output() << "\t" << stPinchSegment_getStart(segment) <<
                    "-" << stPinchSegment_getStart(segment) + 
                    stPinchSegment_getLength(segment) - 1 << " thread #" << 
                    stPinchThread_getName(stPinchSegment_getThread(segment)) <<
                    std::endl;
                
                
            }
            
        }
        
    }
    
    Log::info() << "Counted " << tandemDuplications << " duplications" <<
        std::endl;
    
    // We counted up the tandem duplications. Now return.
    return tandemDuplications;

}

void
writeAdjacencyComponents(
    const std::vector<std::vector<stPinchEnd>>& components,
    const std::string& filename
) {
    // Open file for writing
    std::ofstream out(filename);

    for(auto component: components) {
        out << "Component size " << component.size() << ":" << std::endl;
        
        // For each component, make a set of involved blocks.
        std::unordered_set<stPinchBlock*> blocks;
        
        for(auto end : component) {
            // Get and store the block the end belongs to
            auto block = stPinchEnd_getBlock(&end);
            blocks.insert(block);
            
            // Get the segment the end belongs to
            auto segment = stPinchBlock_getFirst(block);
            
            // Name it. Remember threads are already 1-based.
            std::string segmentName =
                std::to_string(stPinchSegment_getStart(segment)) +
                 "-" + std::to_string(stPinchSegment_getStart(segment) + 
                 stPinchSegment_getLength(segment) - 1);
            
            // And get which end we're on
            char endName = stPinchEnd_getOrientation(&end) ? 'R' : 'L';
            
            // Report the coordinates of that end's block on the reference
            out << "\t" << segmentName << " " << endName << " on thread #" << 
                stPinchThread_getName(stPinchSegment_getThread(segment)) <<
                std::endl;
            
        }
        
        // Do some graphviz
        out << "\tgraph {" << std::endl;
        
        for(auto block : blocks) {
            // Get the segment the end belongs to
            auto segment = stPinchBlock_getFirst(block);
            
            // Go through all the segments in the block in total
            stPinchBlockIt iterator = stPinchBlock_getSegmentIterator(block);
            // And count up the number of them
            size_t count = 0;
            while(stPinchBlockIt_getNext(&iterator) != NULL) {
                count++;
            }
            
            // Make an edge with the info from the first segment, annotated with
            // copy number
            out << "\t\t{rank=same; n" <<
                (stPinchSegment_getStart(segment)) << "L -- n" <<
                std::to_string(stPinchSegment_getStart(segment) + 
                stPinchSegment_getLength(segment)) << 
                "R[color=blue,label=\"" << count << "\"];}" << std::endl;
            
        }
        
        for(auto end : component) {
            // Get the block and segment
            auto block = stPinchEnd_getBlock(&end);
            auto segment = stPinchBlock_getFirst(block);
            
            // Pick left end L or right end R to represent this end
            std::string node = stPinchEnd_getOrientation(&end) ? 
                std::to_string(stPinchSegment_getStart(segment)) + "L" :
                std::to_string(stPinchSegment_getStart(segment) + 
                stPinchSegment_getLength(segment) - 1) + "R";
            
            // Go through all the other ends and connect them to this one
            stSet* otherEnds = stPinchEnd_getConnectedPinchEnds(&end);
            
            // Make an iterator to loop over all the other ends we connect to.
            stSetIterator* iterator = stSet_getIterator(otherEnds); 
            
            stPinchEnd* otherEnd = (stPinchEnd*) stSet_getNext(iterator);
            while(otherEnd != NULL) {
                // For each other end
                
                // Get the block and segment
                auto otherBlock = stPinchEnd_getBlock(otherEnd);
                auto otherSegment = stPinchBlock_getFirst(otherBlock);
                
                // Pick left end L or right end R to represent this other end
                std::string otherNode = stPinchEnd_getOrientation(otherEnd) ? 
                    std::to_string(stPinchSegment_getStart(otherSegment)) +
                    "L" : std::to_string(stPinchSegment_getStart(otherSegment) +
                    stPinchSegment_getLength(otherSegment) - 1) + "R";
                
                if(node < otherNode) {
                    // Only draw these edges going in one direction, so copy
                    // number represents actual copy number.
                    out << "\t\tn" << node << " -- n" << otherNode << 
                        "[color=red];" << std::endl;
                }
                
                otherEnd = (stPinchEnd*) stSet_getNext(iterator);
            }
    
    
            // Clean up sonLib stuff            
            stSet_destructIterator(iterator);
            stSet_destruct(otherEnds);
        }
        
        out << "\t}" << std::endl;
        
    }
    
    
    
    // Close up
    out.close();
}

std::vector<std::vector<stPinchEnd>>
deduplicateIsomorphicAdjacencyComponents(
    const std::vector<std::vector<stPinchEnd>>& components
) {

    // We need these for tracking whether edges are sequence edges or adjacency
    // edges. We put pointers to them on the vflib graph edges.
    static const int SEQUENCE_EDGE = 0;
    static const int ADJACENCY_EDGE = 1;
    
    // We need a vector of vflib Graphs that we can go through in a dumb n^2
    // deduplication way (since we can't really index them any better AFAIK). We
    // also want to keep adjacency components associated with them.
    std::vector<std::pair<Graph*, std::vector<stPinchEnd>>> uniqueSet;
    
    
    // What componsnts will we return?
    std::vector<std::vector<stPinchEnd>> toReturn;
    // Hint to the vector how big it will get, because someone on StackOverflow
    // somewhere said that that was smart.
    // <http://stackoverflow.com/a/7671804/402891>
    toReturn.reserve(uniqueSet.size());
    
    for(auto& pair : uniqueSet) {
        // Go through the final deduplicated set
        
        // Kill all the Graphs left in the deduplicated set.
        delete pair.first;
        
        // Move the component into the new vector, so we can return it.
        toReturn.push_back(std::move(pair.second));
    }
    
    // Return the deduplicated components
    return toReturn;
    
     
    
}


