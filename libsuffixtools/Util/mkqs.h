//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// mkqs - multikey quicksort
//
// Perform a ternary quicksort of strings as described in
// Bentley and Sedgewick, 1997
//
// Example code was downloaded from http://www.cs.princeton.edu/~rs/strings/demo.c
// Modified by JTS to take in a comparator and use a generic type
//

#ifndef MKQS_H
#define MKQS_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <queue>
#include "MkqsThread.h"

#define mkqs_swap(a, b) { T tmp = x[a]; x[a] = x[b]; x[b] = tmp; }

// Swap [i..i+n] and [j..j+n] in x
template<typename T>
inline void vecswap(int i, int j, int n, T* x)
{   
    while (n-- > 0) 
    {
        mkqs_swap(i, j);
        i++;
        j++;
    }
}

#define mkqs_swap2(a, b) { T t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(p) (primarySorter.getChar(*(p), depth))
#define elem2char(e, d) (primarySorter.getChar((e), (d)))

template<typename T>
inline void vecswap2(T* a, T* b, size_t n)
{   
    for(size_t i = 0; i < n; i++)
    {
        T t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

template<typename T, typename PrimarySorter>
T* med3func(T* a, T* b, T* c, int depth, const PrimarySorter& primarySorter)
{   int va, vb, vc;
    if ((va=ptr2char(a)) == (vb=ptr2char(b)))
        return a;
    if ((vc=ptr2char(c)) == va || vc == vb)
        return c;       
    return va < vb ?
          (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}
#define med3(a, b, c) med3func(a, b, c, depth, primarySorter)
template<typename T, typename PrimarySorter, typename FinalSorter>
inline void inssort(T* a, size_t n, size_t d, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{   
    T *pi, *pj, s, t;
    for (pi = a + 1; --n > 0; pi++)
    {
        for (pj = pi; pj > a; pj--) 
        {
            // Inline strcmp: break if *(pj-1) <= *pj
            T elem_s = *(pj - 1);
            T elem_t = *pj;
            const char* s = primarySorter.getChrPtr(elem_s);
            const char* t = primarySorter.getChrPtr(elem_t);

            for (s=s+d, t=t+d; *s==*t && *s!=0; s++, t++)
                ;
            if (*s < *t || (*s == *t && finalSorter(elem_s, elem_t)))
                break;
            mkqs_swap2(pj, pj-1);
        }
    }
}

// Function to audit a block and complain if it isn't sorted.
template<typename T, typename PrimarySorter, typename FinalSorter>
bool checkSort(T* a, size_t n, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{
    for(int64_t i = 0; i < (int64_t)n - 1; i++) {
        // Scan it and assert order.
        
        T& first = a[i];
        T& second = a[i + 1];
        
        // We need a < b

        // First compare the stings.
        int comparison = strcmp(primarySorter.getChrPtr(first),
            primarySorter.getChrPtr(second));
        
        if(comparison > 0) {
            // Complain if the sort by string is wrong.
            std::cerr << "[mkqs] error: incorrect sorting by string!" 
                << std::endl;
        } else if(comparison == 0 && !finalSorter(first, second)) {
            // Complain if the sort by string is right but the sort by index is
            // wrong.
            std::cerr << "[mkqs] error: incorrect sorting by index!" 
                << std::endl;
        } else {
            // Both sorts are correct. Try the next pair of positions.
            continue;
        }
        
        // We already complained that sort is wrong. Say what the sort is
        // actually.
        std::cerr << "Elements:" << std::endl;
        std::cerr << "a[" << i << "] = ID " << first.getID() << 
            " Position " << first.getPos() << std::endl;
        std::cerr << "a[" << i + 1 << "] = ID " << second.getID() << 
            " Position " << second.getPos() << std::endl;
            
        // We found a sort error.
        return false;
            
    }
    
    // No sort errors were found.
    return true;
}

// Serial multi-key quicksort implementation.
template<typename T, typename PrimarySorter, typename FinalSorter>
void mkqs2(T* a, size_t n, size_t depth, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{   
    // r sometimes gets used as a difference between chars, so it needs to be
    // signed. However' it's sometimes used as  a difference between positions,
    // so it needs to be big.
    int64_t r, partval;
    T *pa, *pb, *pc, *pd, *pm, *pn, t;
   
    while(true) {
        // We wrap the entire function in an infinite loop to facilitate manual
        // tail recursion.
   
        if (n < 10) 
        {
            // For small problems, sort with insertion sort. Tell it not to look
            // before the given depth in the strings it's sorting, since we
            // already did our sort up to that depth.
            inssort(a, n, depth, primarySorter, finalSorter);
            return;
        }

        // This line does nothing since we replace pn later.
        pm = a + (n/2);
        // Get a pointer to the last element in our array.
        pn = a + (n-1);

        // Pick a pivot element randomly, and swap it over to the start of the
        // array, wher it will begin a block of elements equal to the pivot.
        size_t mid_idx = rand() % n;
        pm = &a[mid_idx];
        mkqs_swap2(a, pm);

        // Grab the pivot's value.
        partval = ptr2char(a);
        
        // Now we need to partition the array into elements less than the pivot,
        // equal to the pivot, and greater than the pivot.
        
        // Set up places to put things <= than the pivot. pa is at the end of a
        // block of things equal to the pivot, and bp is at the end of a block
        // of things <= the pivot.
        pa = pb = a + 1;
        
        // Set up places to put things >= the pivot. pc is at the start of a
        // block of things >= the pivot, and pd is at the start of a block of
        // things equal to the pivot.
        pc = pd = a + n-1;
        
        // So we're basically splitting up our array like this:
        //   a       pa        pb                   pc          pd     a + n - 1
        //   v       v         v                    v           v          v        
        //  +---------+---------+------------------+-----------+------------+       
        //  | = pivot | < pivot |  unpartitioned   |  > pivot  | = pivot    |       
        //  +---------+---------+------------------+-----------+------------+       
        //                                                                          
        
        
        for (;;) 
        {
            // Sort things into their appropriate blocks, eating away at the
            // unpartitioned space in the middle.
        
            while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) 
            {
                // For things <= the pivot, on the left
            
                if (r == 0) {
                    // Put things equal to the pivot before pa, and advance it
                    // right.
                    mkqs_swap2(pa, pb);
                    pa++;
                }
                
                // Advance pb to encompass one more <= the pivot element.
                pb++;
            }
            while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) 
            {
                // For things >= the pivot, on the right
                
                if (r == 0) {
                    // Put things equal to the pivot at pd, and advance it left.
                    mkqs_swap2(pc, pd);
                    pd--;
                }
                
                // Advance pc to encompass one more >= the pivot element.
                pc--;
            }
            
            // If we've gotten the <= the pivot block on the left to meet the >=
            // the pivot block on the right, we've finished partitioning.
            if (pb > pc) break;
            
            // Otherwise, we got here, because something > the pivot was on the
            // left, and something < the pivot was on the right. so we need to
            // flip those things to be on the correct sides.
            mkqs_swap2(pb, pc);
            // And extend the partition block endpoints inwards accordingly.
            pb++;
            pc--;
        }
        
        // Get a 1-past-the-end pointer.
        pn = a + n;
        
        // Find the smaller of the <-the-pivot and the =-the-pivot blocks on the
        // left, and do the smallest swap needed to bring the =-the-pivot block
        // to the middle.
        r = std::min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
        // Do a similar operation on the right side, so that the order of blocks
        // in the array is now <-the-pivot, =-the-pivot, >-the-pivot.
        r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
        
        // The array is now partitioned.
        
        if ((r = pb-pa) > 1) {
            // There is more than one element less than the pivot, so we need to
            // sort those (on this character). Recurse and do that.
            mkqs2(a, r, depth, primarySorter, finalSorter);
            
            // TODO: If we had a spectacularly bad pivot every time, we could
            // get O(n) recursion depth here.
        }
        if ((r = pd-pc) > 1) {
            // There is more than one element greater than the pivot, so we need
            // to sort those (on this character). Recurse and do that.
            mkqs2(a + n-r, r, depth, primarySorter, finalSorter);
            
            // TODO: If we had a spectacularly bad pivot every time, we could
            // get O(n) recursion depth here.
        }
        
        // How far do we have to go to get to the block of things equal to the
        // pivot?
        r = pb - pa;
        if (ptr2char(a + r) == 0) {
            // The pivot character was '\0', the stop character. The suffixes
            // that were equal to the pivot at this character all end here. So
            // we're done with the primary sort by characters.
            
            // Count up the number of sequences that had the pivot character.
            size_t n2 = pa - a + pn - pd - 1;
            
            // Now that they're sorted by the primary sorter, sort them by the
            // secondary sorter within that.
            std::sort(a + r, a + r + n2, finalSorter);
            
            // Don't do any manually-implemented tail recursion.
            return;
            
        } else {
            // The pivot character was not '\0', the stop character, so the
            // strings equal to the pivot at this character have subsequent
            // characters.
        
            // We need to sort the elements that were equal at this character on
            // the next character, but we need to do it without recursing,
            // because if we recurse at every character we'll fill up the stack
            // when our suffixes are bigger than reads.
            
            // We need to do the equivalent of:
            // mkqs2(a + r, pa-a + pn-pd-1, depth+1, primarySorter,
            //     finalSorter);
            
            // If we could force tail call optimization, we could do that, but
            // we can't. So we loop around.
            
            // Set the arguments
            n = pa - a + pn - pd - 1;
            a = a + r;
            depth++;
            
            // "Recurse" by not returning, and thus looping around.
            
            // TODO: Is there a way to audit the sort here?
        }
    }
}

// Parallel multikey quicksort. It performs mkqs but will
// subdivide the array to sort into sub jobs which can be sorted using threads.
template<typename T, typename PrimarySorter, typename FinalSorter>
void parallel_mkqs(T* pData, size_t n, int numThreads, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{
    typedef MkqsJob<T> Job;
    typedef std::queue<Job> JobQueue;
    Job initialJob(pData, n, 0);
    JobQueue queue;
    queue.push(initialJob);
    
    // Create the mutex that guards the queue
    pthread_mutex_t queue_mutex;
    int ret = pthread_mutex_init(&queue_mutex, NULL);
    if(ret != 0)
    {
        std::cerr << "Mutex initialization failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Calculate the threshold size for performing serial continuation of the sort. Once the chunks 
    // are below this size, it is better to not subdivide the problem into smaller chunks
    // to avoid the overhead of locking, adding to the queue, etc. 
    size_t threshold_size = n / numThreads;

    // Create the semaphore used to signal that data is ready to be processed
    // Initial value is 1 as there is one item on the queue to start
    sem_t queue_sem;
    ret = sem_init( &queue_sem, PTHREAD_PROCESS_PRIVATE, 1 );
    if(ret != 0)
    {
        std::cerr << "Semaphore initialization failed with error " << ret << "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }

    // Create the semaphore indicating a thread is finished working.
    // This semaphore is incremented by the threads in their run loop
    // before waiting for queue items. If the thread takes an item,
    // this semaphore is decremented. The main thread checks
    // this semaphore after confirming the queue is empty. If this 
    // semaphore value equals the total number of threads, 
    // no more work can remain the the threads are cleaned up
    sem_t done_sem;
    ret = sem_init( &done_sem, PTHREAD_PROCESS_PRIVATE, 0 );
    if(ret != 0)
    {
        std::cerr << "Semaphore initialization failed with error " << ret << "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }

    // Create and start the threads
    MkqsThread<T, PrimarySorter, FinalSorter>* threads[numThreads];
    for(int i = 0; i < numThreads; ++i)
    {
        threads[i] = new MkqsThread<T, PrimarySorter, FinalSorter>(i, &queue, &queue_mutex, &queue_sem, &done_sem, threshold_size, &primarySorter, &finalSorter);   
        threads[i]->start();
    }

    // Check for the end condition
    bool done = false;
    while(!done)
    {
        sleep(1);

        // Check if the queue is empty
        // If it is and all threads are finished working (all have posted to
        // the done semaphore), then the threads can be cleaned up.
        pthread_mutex_lock(&queue_mutex);
        if(queue.empty())
        {
            int semval;
            sem_getvalue(&done_sem, &semval);
            if(semval == numThreads)
                done = true;
        }
        pthread_mutex_unlock(&queue_mutex);
    }

    // Signal all the threads to stop, then post to the semaphore they are waiting on
    // All threads will pick up the stop request after the posts and call pthread exit
    for(int i = 0; i < numThreads; ++i)
    {
        threads[i]->stop();
    }

    for(int i = 0; i < numThreads; ++i)
    {
        sem_post(&queue_sem);
    }

    // Join and destroy the threads
    for(int i = 0; i < numThreads; ++i)
    {
        threads[i]->join();
        delete threads[i];
    }

    // Destroy the semaphore
    sem_destroy(&queue_sem);
    sem_destroy(&done_sem);

    // Destroy the queue mutex
    ret = pthread_mutex_destroy(&queue_mutex);
    if(ret != 0)
    {
        std::cerr << "Mutex destruction failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
// Perform a partial sort of the data using the mkqs algorithm
// Iterative sort jobs are created and added to pQueue which is
// protected by pQueueMutex. After addition, pQueueSem is updated.
//
template<typename T, class PrimarySorter, class FinalSorter>
void parallel_mkqs_process(MkqsJob<T>& job, 
                           std::queue<MkqsJob<T> >* pQueue, 
                           pthread_mutex_t* pQueueMutex, 
                           sem_t* pQueueSem,
                           const PrimarySorter& primarySorter, 
                           const FinalSorter& finalSorter)
{
    T* a = job.pData;
    size_t n = job.n;
    size_t depth = job.depth;
    
    int64_t r, partval;
    T *pa, *pb, *pc, *pd, *pm, *pn, t;
    
    if(n < 10) 
    {
        inssort(a, n, depth, primarySorter, finalSorter);
        return;
    }
    
    pm = a + (n/2);
    pn = a + (n-1);

    size_t mid_idx = rand() % n;

    pm = &a[mid_idx];
    mkqs_swap2(a, pm);
    partval = ptr2char(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;
    for (;;) 
    {
        while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) 
        {
            if (r == 0) { mkqs_swap2(pa, pb); pa++; }
            pb++;
        }
        while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) 
        {
            if (r == 0) { mkqs_swap2(pc, pd); pd--; }
            pc--;
        }
        if (pb > pc) break;
        mkqs_swap2(pb, pc);
        pb++;
        pc--;
    }
    pn = a + n;
    r = std::min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);

    // Lock the queue and push new items if necessary
    // If new items are added to the queue, the semaphore is posted to
    pthread_mutex_lock(pQueueMutex);

    if ((r = pb-pa) > 1)
    {
        MkqsJob<T> job(a, r, depth);
        pQueue->push(job);
        sem_post(pQueueSem);
    }
    
    if (ptr2char(a + r) != 0)
    {
        MkqsJob<T> job(a + r, pa-a + pn-pd-1, depth + 1);
        pQueue->push(job);
        sem_post(pQueueSem);
    }
    else
    {
        // Finalize the sort
        size_t n2 = pa - a + pn - pd - 1;
        std::sort(a + r, a + r + n2, finalSorter);
    }

    if ((r = pd-pc) > 1)
    {
        MkqsJob<T> job(a + n-r, r, depth);
        pQueue->push(job);
        sem_post(pQueueSem);
    }

    // Unlock the mutex
    pthread_mutex_unlock(pQueueMutex);
}
#endif

