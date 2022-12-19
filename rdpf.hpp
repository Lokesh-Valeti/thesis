#ifndef __RDPF_HPP__
#define __RDPF_HPP__

#include <vector>
#include <iostream>

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

    RDPF() {}

    // Construct a DPF with the given (XOR-shared) target location, and
    // of the given depth, to be used for random-access memory reads and
    // writes.  The DPF is constructed collaboratively by P0 and P1,
    // with the server P2 helping by providing correlated randomness,
    // such as SelectTriples.
    //
    // Cost:
    // (2 DPFnode + 2 bytes)*depth + 1 word communication in
    // 2*depth + 1 messages
    // (2 DPFnode + 1 byte)*depth communication from P2 to each party
    // 2^{depth+1}-2 local AES operations for P0,P1
    // 0 local AES operations for P2
    RDPF(MPCTIO &tio, yield_t &yield,
        RegXS target, nbits_t depth);

    // The number of bytes it will take to store this RDPF
    size_t size() const;

    // The number of bytes it will take to store a RDPF of the given
    // depth
    static size_t size(nbits_t depth);

    // The depth
    inline nbits_t depth() const { return cw.size(); }
};

// I/O for RDPFs

template <typename T>
T& operator>>(T &is, RDPF &rdpf)
{
    is.read((char *)&rdpf.seed, sizeof(rdpf.seed));
    uint8_t depth;
    assert(depth <= VALUE_BITS);
    is.read((char *)&depth, sizeof(depth));
    rdpf.cw.clear();
    for (uint8_t i=0; i<depth; ++i) {
        DPFnode cw;
        is.read((char *)&cw, sizeof(cw));
        rdpf.cw.push_back(cw);
    }
    value_t cfbits = 0;
    is.read((char *)&cfbits, BITBYTES(depth));
    rdpf.cfbits = cfbits;
    is.read((char *)&rdpf.unit_sum_inverse, sizeof(rdpf.unit_sum_inverse));
    is.read((char *)&rdpf.scaled_sum, sizeof(rdpf.scaled_sum));
    is.read((char *)&rdpf.scaled_xor, sizeof(rdpf.scaled_xor));

    return is;
}

template <typename T>
T& operator<<(T &os, const RDPF &rdpf)
{
    os.write((const char *)&rdpf.seed, sizeof(rdpf.seed));
    uint8_t depth = rdpf.cw.size();
    assert(depth <= VALUE_BITS);
    os.write((const char *)&depth, sizeof(depth));
    for (uint8_t i=0; i<depth; ++i) {
        os.write((const char *)&rdpf.cw[i], sizeof(rdpf.cw[i]));
    }
    os.write((const char *)&rdpf.cfbits, BITBYTES(depth));
    os.write((const char *)&rdpf.unit_sum_inverse, sizeof(rdpf.unit_sum_inverse));
    os.write((const char *)&rdpf.scaled_sum, sizeof(rdpf.scaled_sum));
    os.write((const char *)&rdpf.scaled_xor, sizeof(rdpf.scaled_xor));

    return os;
}

#endif
