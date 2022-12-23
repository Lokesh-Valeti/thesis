#include <bsd/stdlib.h> // arc4random_buf

#include "rdpf.hpp"
#include "bitutils.hpp"
#include "mpcops.hpp"

// Compute the multiplicative inverse of x mod 2^{VALUE_BITS}
// This is the same as computing x to the power of
// 2^{VALUE_BITS-1}-1.
static value_t inverse_value_t(value_t x)
{
    int expon = 1;
    value_t xe = x;
    // Invariant: xe = x^(2^expon - 1) mod 2^{VALUE_BITS}
    // Goal: compute x^(2^{VALUE_BITS-1} - 1)
    while (expon < VALUE_BITS-1) {
        xe = xe * xe * x;
        ++expon;
    }
    return xe;
}

// Construct a DPF with the given (XOR-shared) target location, and
// of the given depth, to be used for random-access memory reads and
// writes.  The DPF is construction collaboratively by P0 and P1,
// with the server P2 helping by providing various kinds of
// correlated randomness, such as MultTriples and AndTriples.
//
// This algorithm is based on Appendix C from the Duoram paper, with a
// small optimization noted below.
RDPF::RDPF(MPCTIO &tio, yield_t &yield,
    RegXS target, nbits_t depth, bool save_expansion)
{
    int player = tio.player();
    size_t &aes_ops = tio.aes_ops();

    // Choose a random seed
    arc4random_buf(&seed, sizeof(seed));
    // Ensure the flag bits (the lsb of each node) are different
    seed = set_lsb(seed, !!player);
    cfbits = 0;
    whichhalf = (player == 1);

    // The root level is just the seed
    nbits_t level = 0;
    DPFnode *curlevel = NULL;
    DPFnode *nextlevel = new DPFnode[1];
    nextlevel[0] = seed;

    // Construct each intermediate level
    while(level < depth) {
        delete[] curlevel;
        curlevel = nextlevel;
        if (save_expansion && level == depth-1) {
            expansion.resize(1<<depth);
            nextlevel = expansion.data();
        } else {
            nextlevel = new DPFnode[1<<(level+1)];
        }
        // Invariant: curlevel has 2^level elements; nextlevel has
        // 2^{level+1} elements

        // The bit-shared choice bit is bit (depth-level-1) of the
        // XOR-shared target index
        RegBS bs_choice = target.bit(depth-level-1);
        size_t curlevel_size = (size_t(1)<<level);
        DPFnode L = _mm_setzero_si128();
        DPFnode R = _mm_setzero_si128();
        // The server doesn't need to do this computation, but it does
        // need to execute mpc_reconstruct_choice so that it sends
        // the AndTriples at the appropriate time.
        if (player < 2) {
            for(size_t i=0;i<curlevel_size;++i) {
                DPFnode lchild, rchild;
                prgboth(lchild, rchild, curlevel[i], aes_ops);
                L = (L ^ lchild);
                R = (R ^ rchild);
                if (nextlevel) {
                    nextlevel[2*i] = lchild;
                    nextlevel[2*i+1] = rchild;
                }
            }
        }
        // If we're going left (bs_choice = 0), we want the correction
        // word to be the XOR of our right side and our peer's right
        // side; if bs_choice = 1, it should be the XOR or our left side
        // and our peer's left side.

        // We also have to ensure that the flag bits (the lsb) of the
        // side that will end up the same be of course the same, but
        // also that the flag bits (the lsb) of the side that will end
        // up different _must_ be different.  That is, it's not enough
        // for the nodes of the child selected by choice to be different
        // as 128-bit values; they also have to be different in their
        // lsb.

        // This is where we make a small optimization over Appendix C of
        // the Duoram paper: instead of keeping separate correction flag
        // bits for the left and right children, we observe that the low
        // bit of the overall correction word effectively serves as one
        // of those bits, so we just need to store one extra bit per
        // level, not two.  (We arbitrarily choose the one for the right
        // child.)

        // Note that the XOR of our left and right child before and
        // after applying the correction word won't change, since the
        // correction word is applied to either both children or
        // neither, depending on the value of the parent's flag. So in
        // particular, the XOR of the flag bits won't change, and if our
        // children's flag's XOR equals our peer's children's flag's
        // XOR, then we won't have different flag bits even for the
        // children that have different 128-bit values.

        // So we compute our_parity = lsb(L^R)^player, and we XOR that
        // into the R value in the correction word computation.  At the
        // same time, we exchange these parity values to compute the
        // combined parity, which we store in the DPF.  Then when the
        // DPF is evaluated, if the parent's flag is set, not only apply
        // the correction work to both children, but also apply the
        // (combined) parity bit to just the right child.  Then for
        // unequal nodes (where the flag bit is different), exactly one
        // of the four children (two for P0 and two for P1) will have
        // the parity bit applied, which will set the XOR of the lsb of
        // those four nodes to just L0^R0^L1^R1^our_parity^peer_parity
        // = 1 because everything cancels out except player (for which
        // one player is 0 and the other is 1).

        bool our_parity_bit = get_lsb(L ^ R) ^ !!player;
        DPFnode our_parity = lsb128_mask[our_parity_bit];

        DPFnode CW;
        bool peer_parity_bit;
        // Exchange the parities and do mpc_reconstruct_choice at the
        // same time (bundled into the same rounds)
        std::vector<coro_t> coroutines;
        coroutines.emplace_back(
            [&](yield_t &yield) {
                tio.queue_peer(&our_parity_bit, 1);
                yield();
                uint8_t peer_parity_byte;
                tio.recv_peer(&peer_parity_byte, 1);
                peer_parity_bit = peer_parity_byte & 1;
            });
        coroutines.emplace_back(
            [&](yield_t &yield) {
                mpc_reconstruct_choice(tio, yield, CW, bs_choice,
                    (R ^ our_parity), L);
            });
        run_coroutines(yield, coroutines);
        bool parity_bit = our_parity_bit ^ peer_parity_bit;
        cfbits |= (value_t(parity_bit)<<level);
        DPFnode CWR = CW ^ lsb128_mask[parity_bit];
        if (player < 2) {
            if (level < depth-1) {
                for(size_t i=0;i<curlevel_size;++i) {
                    bool flag = get_lsb(curlevel[i]);
                    nextlevel[2*i] = xor_if(nextlevel[2*i], CW, flag);
                    nextlevel[2*i+1] = xor_if(nextlevel[2*i+1], CWR, flag);
                }
            } else {
                // Recall there are four potentially useful vectors that
                // can come out of a DPF:
                // - (single-bit) bitwise unit vector
                // - additive-shared unit vector
                // - XOR-shared scaled unit vector
                // - additive-shared scaled unit vector
                //
                // (No single DPF should be used for both of the first
                // two or both of the last two, though, since they're
                // correlated; you _can_ use one of the first two and
                // one of the last two.)
                //
                // For each 128-bit leaf, the low bit is the flag bit,
                // and we're guaranteed that the flag bits (and indeed
                // the whole 128-bit value) for P0 and P1 are the same
                // for every leaf except the target, and that the flag
                // bits definitely differ for the target (and the other
                // 127 bits are independently random on each side).
                //
                // We divide the 128-bit leaf into a low 64-bit word and
                // a high 64-bit word.  We use the low word for the unit
                // vector and the high word for the scaled vector; this
                // choice is not arbitrary: the flag bit in the low word
                // means that the sum of all the low words (with P1's
                // low words negated) across both P0 and P1 is
                // definitely odd, so we can compute that sum's inverse
                // mod 2^64, and store it now during precomputation.  At
                // evaluation time for the additive-shared unit vector,
                // we will output this global inverse times the low word
                // of each leaf, which will make the sum of all of those
                // values 1.  (This technique replaces the protocol in
                // Appendix D of the Duoram paper.)
                //
                // For the scaled vector, we just have to compute shares
                // of what the scaled vector is a sharing _of_, but
                // that's just XORing or adding all of each party's
                // local high words; no communication needed.

                value_t low_sum = 0;
                value_t high_sum = 0;
                value_t high_xor = 0;
                for(size_t i=0;i<curlevel_size;++i) {
                    bool flag = get_lsb(curlevel[i]);
                    DPFnode leftchild = xor_if(nextlevel[2*i], CW, flag);
                    DPFnode rightchild = xor_if(nextlevel[2*i+1], CWR, flag);
                    if (save_expansion) {
                        nextlevel[2*i] = leftchild;
                        nextlevel[2*i+1] = rightchild;
                    }
                    value_t leftlow = value_t(_mm_cvtsi128_si64x(leftchild));
                    value_t rightlow = value_t(_mm_cvtsi128_si64x(rightchild));
                    value_t lefthigh =
                        value_t(_mm_cvtsi128_si64x(_mm_srli_si128(leftchild,8)));
                    value_t righthigh =
                        value_t(_mm_cvtsi128_si64x(_mm_srli_si128(rightchild,8)));
                    low_sum += (leftlow + rightlow);
                    high_sum += (lefthigh + righthigh);
                    high_xor ^= (lefthigh ^ righthigh);
                }
                if (player == 1) {
                    low_sum = -low_sum;
                    high_sum = -high_sum;
                }
                scaled_sum.ashare = high_sum;
                scaled_xor.xshare = high_xor;
                // Exchange low_sum and add them up
                tio.queue_peer(&low_sum, sizeof(low_sum));
                yield();
                value_t peer_low_sum;
                tio.recv_peer(&peer_low_sum, sizeof(peer_low_sum));
                low_sum += peer_low_sum;
                // The low_sum had better be odd
                assert(low_sum & 1);
                unit_sum_inverse = inverse_value_t(low_sum);
            }
            cw.push_back(CW);
        }

        ++level;
    }

    delete[] curlevel;
    if (!save_expansion) {
        delete[] nextlevel;
    }
}

// The number of bytes it will take to store a RDPF of the given depth
size_t RDPF::size(nbits_t depth)
{
    return sizeof(seed) + depth*sizeof(DPFnode) + BITBYTES(depth) +
        sizeof(unit_sum_inverse) + sizeof(scaled_sum) +
        sizeof(scaled_xor);
}

// The number of bytes it will take to store this RDPF
size_t RDPF::size() const
{
    uint8_t depth = cw.size();
    return size(depth);
}

// Get the leaf node for the given input
DPFnode RDPF::leaf(address_t input, size_t &op_counter) const
{
    // If we have a precomputed expansion, just use it
    if (expansion.size()) {
        return expansion[input];
    }

    nbits_t totdepth = depth();
    DPFnode node = seed;
    for (nbits_t d=0;d<totdepth;++d) {
        bit_t dir = !!(input & (address_t(1)<<(totdepth-d-1)));
        node = descend(node, d, dir, op_counter);
    }
    return node;
}

// Expand the DPF if it's not already expanded
//
// This routine is slightly more efficient than repeatedly calling
// Eval::next(), but it uses a lot more memory.
void RDPF::expand(size_t &op_counter)
{
    nbits_t depth = this->depth();
    size_t num_leaves = size_t(1)<<depth;
    if (expansion.size() == num_leaves) return;
    expansion.resize(num_leaves);
    address_t index = 0;
    address_t lastindex = 0;
    DPFnode *path = new DPFnode[depth];
    path[0] = seed;
    for (nbits_t i=1;i<depth;++i) {
        path[i] = descend(path[i-1], i-1, 0, op_counter);
    }
    expansion[index++] = descend(path[depth-1], depth-1, 0, op_counter);
    expansion[index++] = descend(path[depth-1], depth-1, 1, op_counter);
    while(index < num_leaves) {
        // Invariant: lastindex and index will both be even, and
        // index=lastindex+2
        uint64_t index_xor = index ^ lastindex;
        nbits_t how_many_1_bits = __builtin_popcount(index_xor);
        // If lastindex -> index goes for example from (in binary)
        // 010010110 -> 010011000, then index_xor will be
        // 000001110 and how_many_1_bits will be 3.
        // That indicates that path[depth-3] was a left child, and now
        // we need to change it to a right child by descending right
        // from path[depth-4], and then filling the path after that with
        // left children.
        path[depth-how_many_1_bits] =
            descend(path[depth-how_many_1_bits-1],
                depth-how_many_1_bits-1, 1, op_counter);
        for (nbits_t i = depth-how_many_1_bits; i < depth-1; ++i) {
            path[i+1] = descend(path[i], i, 0, op_counter);
        }
        lastindex = index;
        expansion[index++] = descend(path[depth-1], depth-1, 0, op_counter);
        expansion[index++] = descend(path[depth-1], depth-1, 1, op_counter);
    }

    delete[] path;
}

// Construct three RDPFs of the given depth all with the same randomly
// generated target index.
RDPFTriple::RDPFTriple(MPCTIO &tio, yield_t &yield,
    nbits_t depth, bool save_expansion)
{
    // Pick a random XOR share of the target
    xs_target.randomize(depth);

    // Now create three RDPFs with that target, and also convert the XOR
    // shares of the target to additive shares
    std::vector<coro_t> coroutines;
    for (int i=0;i<3;++i) {
        coroutines.emplace_back(
            [&, i](yield_t &yield) {
                dpf[i] = RDPF(tio, yield, xs_target, depth,
                    save_expansion);
            });
    }
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_xs_to_as(tio, yield, as_target, xs_target, depth);
        });
    run_coroutines(yield, coroutines);
}

RDPFTriple::node RDPFTriple::descend(const RDPFTriple::node &parent,
    nbits_t parentdepth, bit_t whichchild,
    size_t &op_counter) const
{
    auto [P0, P1, P2] = parent;
    DPFnode C0, C1, C2;
    C0 = dpf[0].descend(P0, parentdepth, whichchild, op_counter);
    C1 = dpf[1].descend(P1, parentdepth, whichchild, op_counter);
    C2 = dpf[2].descend(P2, parentdepth, whichchild, op_counter);
    return std::make_tuple(C0,C1,C2);
}

RDPFPair::node RDPFPair::descend(const RDPFPair::node &parent,
    nbits_t parentdepth, bit_t whichchild,
    size_t &op_counter) const
{
    auto [P0, P1] = parent;
    DPFnode C0, C1;
    C0 = dpf[0].descend(P0, parentdepth, whichchild, op_counter);
    C1 = dpf[1].descend(P1, parentdepth, whichchild, op_counter);
    return std::make_tuple(C0,C1);
}
