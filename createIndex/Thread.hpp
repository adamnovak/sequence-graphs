#ifndef THREAD_HPP
#define THREAD_HPP

#include <thread>

/**
 * Defines a class that wraps std::thread. Was used to fix up itimers for geprof
 * compatibility, but I can't work out how to forward rvalue references through
 * the crazy bind-thread-lambda stack I built, so it no longer does that.
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
    explicit Thread(_Callable&& __f, _Args&&... __args): 
        std::thread(std::forward<_Callable>(__f), 
        std::forward<_Args>(__args)...) {
        
        // Nothing to do
    }
#endif

};


#endif
