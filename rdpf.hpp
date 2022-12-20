#ifndef __RDPF_HPP__
#define __RDPF_HPP__

#include <vector>
#include <iostream>

#include "mpcio.hpp"
#include "coroutine.hpp"
#include "types.hpp"
#include "bitutils.hpp"

struct RDPF {
    // The 128-bit seed
    DPFnode seed;
    // Which half of the DPF are we?
    bit_t whichhalf;
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
    // If we're saving the expansion, put it here
    std::vector<DPFnode> expansion;

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
        RegXS target, nbits_t depth, bool save_expansion = false);

    // The number of bytes it will take to store this RDPF
    size_t size() const;

    // The number of bytes it will take to store a RDPF of the given
    // depth
    static size_t size(nbits_t depth);

    // The depth
    inline nbits_t depth() const { return cw.size(); }

    // Descend from a node at depth parentdepth to one of its children
    // whichchild = 0: left child
    // whichchild = 1: right child
    //
    // Cost: 1 AES operation
    DPFnode descend(const DPFnode parent, nbits_t parentdepth,
        bit_t whichchild, size_t &op_counter) const;

    // Get the leaf node for the given input
    //
    // Cost: depth AES operations
    DPFnode leaf(address_t input, size_t &op_counter) const;

    // Expand the DPF if it's not already expanded
    void expand(size_t &op_counter);

    // Streaming evaluation, to avoid taking up enough memory to store
    // an entire evaluation
    class Eval {
        friend class RDPF; // So eval() can call the Eval constructor
        const RDPF &rdpf;
        size_t &op_counter;
        bool use_expansion;
        nbits_t depth;
        address_t indexmask;
        address_t pathindex;
        address_t nextindex;
        std::vector<DPFnode> path;
        Eval(const RDPF &rdpf, size_t &op_counter, address_t start,
            bool use_expansion);
    public:
        DPFnode next();
    };

    // Create an Eval object that will start its output at index start.
    // It will wrap around to 0 when it hits 2^depth.  If use_expansion
    // is true, then if the DPF has been expanded, just output values
    // from that.  If use_expansion=false or if the DPF has not been
    // expanded, compute the values on the fly.
    Eval eval(address_t start, size_t &op_counter,
        bool use_expansion=true) const;

    // Get the bit-shared unit vector entry from the leaf node
    inline RegBS unit_bs(DPFnode leaf) const {
        RegBS b;
        b.bshare = get_lsb(leaf);
        return b;
    }

    // Get the additive-shared unit vector entry from the leaf node
    inline RegAS unit_as(DPFnode leaf) const {
        RegAS a;
        value_t lowword = value_t(_mm_cvtsi128_si64x(leaf));
        if (whichhalf == 1) {
            lowword = -lowword;
        }
        a.ashare = lowword * unit_sum_inverse;
        return a;
    }

    // Get the XOR-shared scaled vector entry from the leaf ndoe
    inline RegXS scaled_xs(DPFnode leaf) const {
        RegXS x;
        value_t highword =
            value_t(_mm_cvtsi128_si64x(_mm_srli_si128(leaf,8)));
        x.xshare = highword;
        return x;
    }

    // Get the additive-shared scaled vector entry from the leaf ndoe
    inline RegAS scaled_as(DPFnode leaf) const {
        RegAS a;
        value_t highword =
            value_t(_mm_cvtsi128_si64x(_mm_srli_si128(leaf,8)));
        if (whichhalf == 1) {
            highword = -highword;
        }
        a.ashare = highword;
        return a;
    }

};

// I/O for RDPFs

template <typename T>
T& operator>>(T &is, RDPF &rdpf)
{
    is.read((char *)&rdpf.seed, sizeof(rdpf.seed));
    uint8_t depth;
    // The whichhalf bit is the high bit of depth
    is.read((char *)&depth, sizeof(depth));
    rdpf.whichhalf = !!(depth & 0x80);
    depth &= 0x7f;
    bool read_expanded = false;
    if (depth > 64) {
        read_expanded = true;
        depth -= 64;
    }
    assert(depth <= ADDRESS_MAX_BITS);
    rdpf.cw.clear();
    for (uint8_t i=0; i<depth; ++i) {
        DPFnode cw;
        is.read((char *)&cw, sizeof(cw));
        rdpf.cw.push_back(cw);
    }
    if (read_expanded) {
        rdpf.expansion.resize(1<<depth);
        is.read((char *)rdpf.expansion.data(),
            sizeof(rdpf.expansion[0])<<depth);
    }
    value_t cfbits = 0;
    is.read((char *)&cfbits, BITBYTES(depth));
    rdpf.cfbits = cfbits;
    is.read((char *)&rdpf.unit_sum_inverse, sizeof(rdpf.unit_sum_inverse));
    is.read((char *)&rdpf.scaled_sum, sizeof(rdpf.scaled_sum));
    is.read((char *)&rdpf.scaled_xor, sizeof(rdpf.scaled_xor));

    return is;
}

// Write the DPF to the output stream.  If expanded=true, then include
// the expansion _if_ the DPF is itself already expanded.  You can use
// this to write DPFs to files.
template <typename T>
T& write_maybe_expanded(T &os, const RDPF &rdpf,
    bool expanded = true)
{
    os.write((const char *)&rdpf.seed, sizeof(rdpf.seed));
    uint8_t depth = rdpf.cw.size();
    assert(depth <= ADDRESS_MAX_BITS);
    // The whichhalf bit is the high bit of depth
    // If we're writing an expansion, add 64 to depth as well
    uint8_t whichhalf_and_depth = depth |
        (uint8_t(rdpf.whichhalf)<<7);
    bool write_expansion = false;
    if (expanded && rdpf.expansion.size() == (size_t(1)<<depth)) {
        write_expansion = true;
        whichhalf_and_depth += 64;
    }
    os.write((const char *)&whichhalf_and_depth,
        sizeof(whichhalf_and_depth));
    for (uint8_t i=0; i<depth; ++i) {
        os.write((const char *)&rdpf.cw[i], sizeof(rdpf.cw[i]));
    }
    if (write_expansion) {
        os.write((const char *)rdpf.expansion.data(),
            sizeof(rdpf.expansion[0])<<depth);
    }
    os.write((const char *)&rdpf.cfbits, BITBYTES(depth));
    os.write((const char *)&rdpf.unit_sum_inverse, sizeof(rdpf.unit_sum_inverse));
    os.write((const char *)&rdpf.scaled_sum, sizeof(rdpf.scaled_sum));
    os.write((const char *)&rdpf.scaled_xor, sizeof(rdpf.scaled_xor));

    return os;
}

// The ordinary << version never writes the expansion, since this is
// what we use to send DPFs over the network.
template <typename T>
T& operator<<(T &os, const RDPF &rdpf)
{
    return write_maybe_expanded(os, rdpf, false);
}

// Computational peers will generate triples of RDPFs with the _same_
// random target for use in Duoram.  They will each hold a share of the
// target (neither knowing the complete target index).  They will each
// give one of the DPFs (not a matching pair) to the server, but not the
// shares of the target index.  So computational peers will hold a
// RDPFTriple (which includes both an additive and an XOR share of the
// target index), while the server will hold a RDPFPair (which does
// not).

struct RDPFTriple {
    RegAS as_target;
    RegXS xs_target;
    RDPF dpf[3];

    RDPFTriple() {}

    // Construct three RDPFs of the given depth all with the same
    // randomly generated target index.
    RDPFTriple(MPCTIO &tio, yield_t &yield,
        nbits_t depth, bool save_expansion = false);
};

// I/O for RDPF Triples

// We never write RDPFTriples over the network, so always write
// the DPF expansions if they're available.
template <typename T>
T& operator<<(T &os, const RDPFTriple &rdpftrip)
{
    write_maybe_expanded(os, rdpftrip.dpf[0], true);
    write_maybe_expanded(os, rdpftrip.dpf[1], true);
    write_maybe_expanded(os, rdpftrip.dpf[2], true);
    nbits_t depth = rdpftrip.dpf[0].depth();
    os.write((const char *)&rdpftrip.as_target.ashare, BITBYTES(depth));
    os.write((const char *)&rdpftrip.xs_target.xshare, BITBYTES(depth));
    return os;
}

template <typename T>
T& operator>>(T &is, RDPFTriple &rdpftrip)
{
    is >> rdpftrip.dpf[0] >> rdpftrip.dpf[1] >> rdpftrip.dpf[2];
    nbits_t depth = rdpftrip.dpf[0].depth();
    rdpftrip.as_target.ashare = 0;
    is.read((char *)&rdpftrip.as_target.ashare, BITBYTES(depth));
    rdpftrip.xs_target.xshare = 0;
    is.read((char *)&rdpftrip.xs_target.xshare, BITBYTES(depth));
    return is;
}

struct RDPFPair {
    RDPF dpf[2];
};

// I/O for RDPF Pairs

// We never write RDPFPairs over the network, so always write
// the DPF expansions if they're available.
template <typename T>
T& operator<<(T &os, const RDPFPair &rdpfpair)
{
    write_maybe_expanded(os, rdpfpair.dpf[0], true);
    write_maybe_expanded(os, rdpfpair.dpf[1], true);
    return os;
}

template <typename T>
T& operator>>(T &is, RDPFPair &rdpfpair)
{
    is >> rdpfpair.dpf[0] >> rdpfpair.dpf[1];
    return is;
}

#endif
