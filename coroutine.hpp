#ifndef __COROUTINE_HPP__
#define __COROUTINE_HPP__

#include <vector>
#include <boost/coroutine2/coroutine.hpp>

typedef boost::coroutines2::coroutine<void>::pull_type  coro_t;
typedef boost::coroutines2::coroutine<void>::push_type  yield_t;

// Use this version from the top level; i.e., not running in a coroutine
// itself
inline void run_coroutines(MPCTIO &tio, std::vector<coro_t> &coroutines) {
    // Loop until all the coroutines are finished
    bool finished = false;
    while(!finished) {
        // If this current function is not itself a coroutine (i.e.,
        // this is the top-level function that launches all the
        // coroutines), here's where to call send().  Otherwise, call
        // yield() here to let other coroutines at this level run.
        tio.send();
        finished = true;
        for (auto &c : coroutines) {
            // This tests if coroutine c still has work to do (is not
            // finished)
            if (c) {
                finished = false;
                // Resume coroutine c from the point it yield()ed
                c();
            }
        }
    }
}

// Use this version if a coroutine needs to itself use multiple
// subcoroutines.
inline void run_coroutines(yield_t &yield, std::vector<coro_t> &coroutines) {
    // Loop until all the coroutines are finished
    bool finished = false;
    while(!finished) {
        // If this current function is not itself a coroutine (i.e.,
        // this is the top-level function that launches all the
        // coroutines), here's where to call send().  Otherwise, call
        // yield() here to let other coroutines at this level run.
        yield();
        finished = true;
        for (auto &c : coroutines) {
            // This tests if coroutine c still has work to do (is not
            // finished)
            if (c) {
                finished = false;
                // Resume coroutine c from the point it yield()ed
                c();
            }
        }
    }
}

#endif
