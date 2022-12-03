#ifndef __COROUTINE_HPP__
#define __COROUTINE_HPP__

#include <vector>
#include <boost/coroutine2/coroutine.hpp>

typedef boost::coroutines2::coroutine<void>::pull_type  coro_t;
typedef boost::coroutines2::coroutine<void>::push_type  yield_t;

inline void run_coroutines(MPCIO &mpcio, std::vector<coro_t> &coroutines) {
    // Loop until all the coroutines are finished
    bool finished = false;
    while(!finished) {
        // If this current function is not itself a coroutine (i.e.,
        // this is the top-level function that launches all the
        // coroutines), here's where to call send().  Otherwise, call
        // yield() here to let other coroutines at this level run.
        mpcio.sendall();
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
