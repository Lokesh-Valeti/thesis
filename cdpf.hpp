#ifndef __CDPF_HPP__
#define __CDPF_HPP__

#include <tuple>
#include "types.hpp"
#include "dpf.hpp"

// DPFs for doing comparisons of (typically) 64-bit values. We use the
// technique from:
//
// Kyle Storrier, Adithya Vadapalli, Allan Lyons, Ryan Henry.
// Grotto: Screaming fast (2 + 1)-PC for Z_{2^n} via (2, 2)-DPFs
//
// The idea is that we have a pair of DPFs with 64-bit inputs and a
// single-bit output.  The outputs of these DPFs are the same for all
// 64-bit inputs x except for one special one (target), where they're
// different, but if you have just one of the DPFs, you can't tell what
// the value of target is.  The construction of the DPF is a binary
// tree, where each interior node has a 128-bit value, the low bit of
// which is the "flag" bit.  The invariant is that if a node is on the
// path leading to the target, then not only are the two 128-bit values
// on the node (one from each DPF) different, but their flag (low) bits
// are themselves different, and if a node is not on the path leading to
// the target, then its 128-bit value is the _same_ in the two DPFs.
// Each DPF also comes with an additive share (target0 or target1) of
// the random target value.
//
// Given additive shares x0 and x1 of x, two parties can determine
// bitwise shares of whether x>0 as follows: exchange (target0-x0) and
// (target1-x1); both sides add them to produce S = (target-x).
// Notionally consider (but do not actually construct) a bit vector V of
// length 2^64 with 1s at positions S+1, S+2, ..., S+(2^63-1), wrapping
// around if the indices exceed 2^64-1.  Now consider (but again do not
// actually do) the dot product of V with the full evaluation of the
// DPFs.  The full evaluations of the DPFs are random bit vectors that
// differ in only the bit at position target, so the two dot products
// (which are each a single bit) will be a bitwise shraring of the value
// of V at position target.  Note that if V[target] = 1, then target =
// S+k for some 1 <= k <= 2^63-1, then since target = S+x, we have that
// x = k is in that same range; i.e. x>0 as a 64-bit signed integer (and
// similarly if V[target] = 0, then x <= 0.
//
// So far, this is all standard, and for DPFs of smaller depth, this is
// the same technique we're doing for RDPFs.  But we can't do it for
// vectors of size 2^64; that's too big.  Even for 2^32 it would be
// annoying.  The observation made in the Grotto paper is that you can
// actually compute this bit sharing in time linear in the *depth* of
// the DPF (that is, logarithmic in the length of V), for some kinds of
// vectors V, including the "single block of 1s" one described above.
//
// The key insight is that if you look at any _interior_ node of the
// tree, the corresponding nodes on the two DPFs will be a bit sharing
// of the sum of all the leaves in the subtree rooted at that interior
// node: 0 if target is not in that subtree, and 1 if it is.  So you
// just have to find the minimal set of interior nodes such that the
// leaves of the subtrees rooted at those nodes is exactly the block of
// 1s in V, and then each party adds up the flag bits of those leaves.
// The result is a bit sharing of 1 if V[target]=1 and 0 if V[target]=0;
// that is, it is a bit sharing of V[target], and so (as above) of the
// result of the comparison [x>0].  You can also find and evaluate the
// flag bits of this minimal set in time and memory linear in the depth
// of the DPF.
//
// So at the end, we've computed a bit sharing of [x>0] with local
// computation linear in the depth of the DPF (concretely, fewer than
// 200 AES operations), and only a *single word* of communication in
// each direction (exchanging the target{i}-x{i} values).  Of course,
// this assumes you have one pair of these DPFs lying around, and you
// have to use a fresh pair with a fresh random target value for each
// comparison, since revealing target-x for two different x's but the
// same target leaks the difference of the x's. But in the 3-party
// setting (or even the 2+1-party setting), you can just have the server
// precompute a bunch of these pairs in advance, and hand bunches of the
// first item in each pair to player 0 and the second item in each pair
// to player 1, at preprocessing time (a single message from the server
// to each of player 0 and player 1), and these DPFs are very fast to
// compute, and very small (< 1KB each) to transmit and store.

// See also dpf.hpp for the differences between these DPFs and the ones
// we use for oblivious random access to memory.

struct CDPF : public DPF {
    // Additive and XOR shares of the target value
    RegAS as_target;
    RegXS xs_target;
    // The extra correction word we'll need for the right child at the
    // final leaf layer; this is needed because we're making the tree 7
    // layers shorter than you would naively expect (depth 57 instead of
    // 64), and having the 128-bit labels on the leaf nodes directly
    // represent the 128 bits that would have come out of the subtree of
    // a (notional) depth-64 tree rooted at that depth-57 node.
    DPFnode leaf_cwr;

    // Generate a pair of CDPFs with the given target value
    static std::tuple<CDPF,CDPF> generate(value_t target, size_t &aes_ops);

    // Generate a pair of CDPFs with a random target value
    static std::tuple<CDPF,CDPF> generate(size_t &aes_ops);

    // Descend from the parent of a leaf node to the leaf node
    inline DPFnode descend_to_leaf(const DPFnode &parent,
        bit_t whichchild, size_t &aes_ops) const;
};

// Descend from the parent of a leaf node to the leaf node
inline DPFnode CDPF::descend_to_leaf(const DPFnode &parent,
    bit_t whichchild, size_t &aes_ops) const
{
    DPFnode prgout;
    bool flag = get_lsb(parent);
    prg(prgout, parent, whichchild, aes_ops);
    if (flag) {
        DPFnode CW = cw.back();
        DPFnode CWR = leaf_cwr;
        prgout ^= (whichchild ? CWR : CW);
    }
    return prgout;
}

#include "cdpf.tcc"

#endif
