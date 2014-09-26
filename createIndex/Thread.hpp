#ifndef THREAD_HPP
#define THREAD_HPP

#include <thread>
#include <sys/time.h>
#include <errno.h>

/**
 * Defines a class that wraps std::thread and preserves the itimer of the
 * original thread. Necessary for gprof to profile a multithreaded application.
 * Also provides a convenient hook to make the program single-threaded.
 */
class Thread : public std::thread {

public:
   
#ifdef SINGLE_THREAD
    /**
     * This mirrors the std::thread constructor exactly, but runs the thread
     * contents right now.
     */
    template<typename _Callable, typename... _Args>
    explicit Thread(_Callable&& __f, _Args&&... __args): std::thread([]() {}) {
        
        // Run right away.    
        std::bind(__f, __args...)();
        
    }
#else
    /**
     * This mirrors the std::thread constructor exactly.
     */
    template<typename _Callable, typename... _Args>
    explicit Thread(_Callable&& __f, 
        _Args&&... __args): 
        // We become a new thread
        std::thread(
            // Which executes a closure
            std::bind(
                // Which when executed runs this lambda in the new thread
                [](struct itimerval itimer, _Callable& f, _Args&... args) {
            
                    // Which first sets up the timer.
                    if(setitimer(ITIMER_PROF, &itimer, NULL)) {
                        perror("Can't set itimer!");
                        throw std::runtime_error("Can't set itimer!");
                    }
                    
                    // And then runs the callable on the args.
                    std::bind(f, args...)();
                    
        
            
                },
                // With these arguments (now lvales inside the closure which is
                // itself moved to the thread)
                getitimer_adapter(), std::forward<_Callable>(__f), 
                std::forward<_Args>(__args)...
            )
        ) {
        // Nothing to do here. Already initialized the thread with a lambda in
        // the initializer list.
    }
#endif
    
private:

    /**
     * Private function to grab the itimer value.
     */
    inline struct itimerval getitimer_adapter() {
        // have a local
        struct itimerval value;
        
        // Grab the value by pointer.
        if(getitimer(ITIMER_PROF, &value)) {
            perror("Can't get itimer!");
            throw std::runtime_error("Can't get itimer!");
        }
        
        // Return it by value
        return value;
    }


};


#endif
