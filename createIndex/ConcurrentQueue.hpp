#ifndef CONCURRENTQUEUE_HPP
#define CONCURRENTQUEUE_HPP

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

/**
 * A queue which comes with a lock for controlling access from multiple threads.
 * C++11 only.
 * 
 * See <http://juanchopanzacpp.wordpress.com/2013/02/26/concurrent-queue-c11/>
 */
 
template <typename T>
class ConcurrentQueue {

public:
    
    /**
     * This defines a lock we can give or demand for access.
     */
    typedef std::unique_lock<std::mutex> Lock;
    
    /**
     * Construct a ConcurrentQueue that doesn't track writers (i.e. a normal
     * one).
     */
    ConcurrentQueue(): queue(), mutex(), nonempty(), numWriters() {
    }
    
    /**
     * Construct a ConcurrentQueue that tracks writers. Takes a number of
     * writers, and expects each of those writers to eventually call close().
     * When all the writers have called close(), readers will be able to find
     * out and not block on data from a queue that nobody will ever write to.
     */
    ConcurrentQueue(size_t numWriters): queue(), mutex(), nonempty(), 
        numWriters(numWriters) {
    
    }
    
    /**
     * Lock the queue so it can be safely read or written to. Caller must not
     * already hold a lock on the queue.
     */
    Lock lock() {
        return std::unique_lock<std::mutex>(mutex);
    }
    
    /**
     * Remove and return the first thing in the queue. Caller must hold a lock
     * on the queue, and the queue must be nonempty. The lock passed is
     * released.
     *
     * TODO: would it be faster to hold onto the lock and check again for stuff?
     */
    T dequeue(Lock& callerLock) {
        // Grab the first element.
        T toReturn = queue.front();
        
        // Remove it form our queue.
        queue.pop();
        
        // Unlock the caller's lock.
        callerLock.unlock();
        
        // Give the value to the caller.
        return toReturn;
    }
    
    /**
     * Add something to the end of the queue. Caller must hold a lock on the
     * queue.
     *
     * The lock passed is released, and a waiting thread, if any, is notified
     * that data is available.
     */
    void enqueue(const T& value, Lock& callerLock) {
        // Put the element at the end of the queue.
        queue.push(value);
        
        // Release the caller's lock
        callerLock.unlock();
        
        // Tell anyone who is waiting.
        nonempty.notify_one();
    }
    
    /**
     * Returns true if the queue is empty, and false otherwise. Caller must hold
     * a lock on the queue, but it is not released.
     */
    bool isEmpty(Lock& callerLock) {
        return queue.empty();
    }
    
    /**
     * Wait for the queue to be nonempty and lock it. Will return a lock on a
     * nonempty queue. Caller must not already hold a lock on the queue.
     */
    Lock waitForNonempty() {
        // Lock the queue
        Lock waitLock = lock();
        while(isEmpty()) {
            // If the queue is empty, unlock it and wait for someone to ring the
            // nonempty bell. When that happens, re-lock it and check again to
            // see if we won the item.
            nonempty.wait(waitLock);
        }
        // When we get here, we successfully grabbed the queue with a thing in
        // it.
        return waitLock;
    }
    
    /**
     * Close the queue. Should be called from a writer that the queue was told
     * about in the constructor. That writer must hold a lock on the queue. Each
     * writer must call close() exactly once if the writer counting feature is
     * being used at all. After calling this function, a thread my not write to
     * the queue anymore.
     *
     * The lock passed is released.
     */
    void close(Lock& callerLock) {
        // Say we're done.
        numWriters--;
        
        // Release the lock
        callerLock.unlock();
        
        if(numWriters == 0) {
            // We were the last writer. Tell all the things waiting for data to
            // check again
            nonempty.notify_all();
        }
        
        Log::info() << "Queue writers remaining: " << numWriters << std::endl;
    }
    
    /**
     * Wait for either the queue to become nonempty, or for there to be no more
     * open writers left (i.e. no possibility of any data ever coming). Caller
     * must not already hold a lock on the queue. Only works if the queue was
     * told that a nonzero number of writers exist in the constructor.
     */
    Lock waitForNonemptyOrEnd() {
        // Lock the queue
        Lock waitLock = lock();
        while(isEmpty(waitLock) && numWriters > 0) {
            // If the queue is empty, unlock it and wait for someone to ring the
            // nonempty bell. When that happens, re-lock it and check again to
            // see if we won the item.
            nonempty.wait(waitLock);
        }
        // When we get here, we successfully grabbed the queue with a thing in
        // it.
        return waitLock;
    }
    
    /**
     * Ask for the default move assignment operator.
     */
    ConcurrentQueue<T>& operator=(ConcurrentQueue<T>&& other) = default;
    
    /**
     * Ask for the default move constructor.
     */
    ConcurrentQueue(ConcurrentQueue<T>&& other) = default;

protected:
    // Keep around an actual queue that we controll access to.
    std::queue<T> queue;
    
    // This mutex controls access to the queue, which only one thread can have
    // at a time.
    std::mutex mutex;
    
    // This condition variable lets us signal that the queue is probably non-
    // empty, if someone was waiting for that.
    std::condition_variable nonempty;
    
    // We would have a condition for non-full, but instead we let our queue grow
    // indefinitely. Don't fill up memory.
    
    // How many writers are writing to the queue still (if writer counting is
    // enabled).
    size_t numWriters;
    
private:
    
    /**
     * No copy constructor is allowed.
     */
    ConcurrentQueue<T>(const ConcurrentQueue<T>& other) = delete;
    
    /**
     * No assignment operator either.
     */
    ConcurrentQueue<T>& operator=(const ConcurrentQueue<T> other) = delete;
    
    // Move assignment and move are OK
     
    

};

#endif
