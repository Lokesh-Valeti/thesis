#ifndef __RDPF_HPP__
#define __RDPF_HPP__

#include <vector>

#include "mpcio.hpp"
#include "coroutine.hpp"
#include "types.hpp"

struct RDPF {
    // The 128-bit seed
    DPFnode seed;
    // correction words; the depth of the DPF is the length of this
    // vector
    std::vector<DPFnode> cw;
    // correction flag bits: the one for level i is bit i of this word
    value_t cfbits;

    // Construct a DPF with the given (XOR-shared) target location, and
    // of the given depth, to be used for random-access memory reads and
    // writes.  The DPF is construction collaboratively by P0 and P1,
    // with the server P2 helping by providing various kinds of
    // correlated randomness, such as MultTriples and AndTriples.
    RDPF(MPCTIO &tio, yield_t &yield,
        RegXS target, nbits_t depth);
};

#endif
