#ifndef __RDPF_HPP__
#define __RDPF_HPP__

#include "mpcio.hpp"
#include "coroutine.hpp"
#include "types.hpp"

// Construct a DPF of the given depth to be used for random-access
// memory reads and writes.  The DPF is construction collaboratively by
// P0 and P1, with the server P2 helping by providing various kinds of
// correlated randomness, such as MultTriples and AndTriples.
void rdpf_gen(MPCTIO &tio, yield_t &yield,
    RDPF &rdpf, nbits_t depth);

#endif
