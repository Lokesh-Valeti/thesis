#ifndef __COROUTINE_HPP__
#define __COROUTINE_HPP__

#include <vector>
#include <boost/coroutine2/coroutine.hpp>

typedef boost::coroutines2::coroutine<void>::pull_type  coro_t;
typedef boost::coroutines2::coroutine<void>::push_type  yield_t;

// The top-level coroutine runner will call run_coroutines with
// a MPCTIO, and we should call its send() method.  Subcoroutines that
// launch their own coroutines (and coroutine running) will call
// run_coroutines with a yield_t instead, which we should just call, in
// order to yield to the next higher level of coroutine runner.
static inline void send_or_yield(MPCTIO &tio) { tio.send(); }
static inline void send_or_yield(yield_t &yield) { yield(); }

template <typename T>
inline void run_coroutines(T &mpctio_or_yield, std::vector<coro_t> &coroutines) {
    // Loop until all the coroutines are finished
    bool finished = false;
    while(!finished) {
        // If this current function is not itself a coroutine (i.e.,
        // this is the top-level function that launches all the
        // coroutines), here's where to call send().  Otherwise, call
        // yield() here to let other coroutines at this level run.
        send_or_yield(mpctio_or_yield);
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
