#include <bsd/stdlib.h> // arc4random_buf

#include "rdpf.hpp"
#include "bitutils.hpp"
#include "mpcops.hpp"
#include "aes.hpp"
#include "prg.hpp"

#ifdef DPF_DEBUG
static void dump_node(DPFnode node, const char *label = NULL)
{
    if (label) printf("%s: ", label);
    for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&node)[15-i]); } printf("\n");
}

static void dump_level(DPFnode *nodes, size_t num, const char *label = NULL)
{
    if (label) printf("%s:\n", label);
    for (size_t i=0;i<num;++i) {
        dump_node(nodes[i]);
    }
    printf("\n");
}
#endif

// Construct a DPF with the given (XOR-shared) target location, and
// of the given depth, to be used for random-access memory reads and
// writes.  The DPF is construction collaboratively by P0 and P1,
// with the server P2 helping by providing various kinds of
// correlated randomness, such as MultTriples and AndTriples.
//
// This algorithm is based on Appendix C from the Duoram paper, with a
// small optimization noted below.
RDPF::RDPF(MPCTIO &tio, yield_t &yield,
    RegXS target, nbits_t depth)
{
    int player = tio.player();
    size_t &aesops = tio.aes_ops();

    // Choose a random seed
    arc4random_buf(&seed, sizeof(seed));
    // Ensure the flag bits (the lsb of each node) are different
    seed = set_lsb(seed, !!player);
    cfbits = 0;

    // The root level is just the seed
    nbits_t level = 0;
    DPFnode *curlevel = NULL;
    DPFnode *nextlevel = new DPFnode[1];
    nextlevel[0] = seed;

    // Construct each intermediate level
    while(level < depth) {
        delete[] curlevel;
        curlevel = nextlevel;
        // We don't need to store the last level
        if (level < depth-1) {
            nextlevel = new DPFnode[1<<(level+1)];
        } else {
            nextlevel = NULL;
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
                prgboth(lchild, rchild, curlevel[i], aesops);
                L = _mm_xor_si128(L, lchild);
                R = _mm_xor_si128(R, rchild);
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

        // We have to ensure that the flag bits (the lsb) of the side
        // that will end up the same be of course the same, but also
        // that the flag bits (the lsb) of the side that will end up
        // different _must_ be different.  That is, it's not enough for
        // the nodes of the child selected by choice to be different as
        // 128-bit values; they also have to be different in their lsb.

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

        bool our_parity_bit = get_lsb(_mm_xor_si128(L,R)) ^ !!player;
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
                tio.recv_peer(&peer_parity_bit, 1);
            });
        coroutines.emplace_back(
            [&](yield_t &yield) {
                mpc_reconstruct_choice(tio, yield, CW, bs_choice,
                    _mm_xor_si128(R,our_parity), L);
            });
        run_coroutines(yield, coroutines);
        bool parity_bit = our_parity_bit ^ peer_parity_bit;
        cfbits |= (size_t(parity_bit)<<level);
        DPFnode CWR = _mm_xor_si128(CW,lsb128_mask[parity_bit]);
        if (player < 2) {
            if (nextlevel) {
                for(size_t i=0;i<curlevel_size;++i) {
                    bool flag = get_lsb(curlevel[i]);
                    nextlevel[2*i] = xor_if(nextlevel[2*i], CW, flag);
                    nextlevel[2*i+1] = xor_if(nextlevel[2*i+1], CWR, flag);
                }
            }
            cw.push_back(CW);
        }

        ++level;
    }

    delete[] curlevel;
    delete[] nextlevel;
}
