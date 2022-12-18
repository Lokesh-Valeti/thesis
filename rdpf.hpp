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
    // The amount we have to scale the low words of the leaf values by
    // to get additive shares of a unit vector
    value_t unit_sum_inverse;
    // Additive share of the scaling value M_as such that the high words
    // of the leaf values for P0 and P1 add to M_as * e_{target}
    RegAS scaled_sum;
    // XOR share of the scaling value M_xs such that the high words
    // of the leaf values for P0 and P1 XOR to M_xs * e_{target}
    RegXS scaled_xor;

    // Construct a DPF with the given (XOR-shared) target location, and
    // of the given depth, to be used for random-access memory reads and
    // writes.  The DPF is construction collaboratively by P0 and P1,
    // with the server P2 helping by providing various kinds of
    // correlated randomness, such as MultTriples and AndTriples.
    //
    // Cost:
    // (3 DPFnode + 1 byte)*depth + 1 word communication in
    // 2*depth + 1 messages
    // 3*depth DPFnode communication from P2 to each party
    // 2^{depth+1}-2 local AES operations
    RDPF(MPCTIO &tio, yield_t &yield,
        RegXS target, nbits_t depth);
};

#endif
