#include <bsd/stdlib.h> // arc4random_buf

#include "rdpf.hpp"
#include "bitutils.hpp"
#include "mpcops.hpp"
#include "aes.hpp"
#include "prg.hpp"

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

// Construct a DPF with the given (XOR-shared) target location, and
// of the given depth, to be used for random-access memory reads and
// writes.  The DPF is construction collaboratively by P0 and P1,
// with the server P2 helping by providing various kinds of
// correlated randomness, such as MultTriples and AndTriples.
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
    while(level < depth - 1) {
        delete[] curlevel;
        curlevel = nextlevel;
        nextlevel = new DPFnode[1<<(level+1)];
        // Invariant: curlevel has 2^level elements; nextlevel has
        // 2^{level+1} elements

        // The bit-shared choice bit is bit (depth-level-1) of the
        // XOR-shared target index
        RegBS bs_choice = target.bit(depth-level-1);
        size_t curlevel_size = (size_t(1)<<level);
        DPFnode L = _mm_setzero_si128();
        DPFnode R = _mm_setzero_si128();
        if (player < 2) {
            for(size_t i=0;i<curlevel_size;++i) {
                prgboth(nextlevel[2*i], nextlevel[2*i+1], curlevel[i], aesops);
                L = _mm_xor_si128(L, nextlevel[2*i]);
                R = _mm_xor_si128(R, nextlevel[2*i+1]);
            }
        }
        DPFnode CW;
        mpc_reconstruct_choice(tio, yield, CW, bs_choice, R, L);
        if (player < 2) {
            for(size_t i=0;i<curlevel_size;++i) {
                bool flag = get_lsb(curlevel[i]);
                nextlevel[2*i] = xor_if(nextlevel[2*i], CW, flag);
                nextlevel[2*i+1] = xor_if(nextlevel[2*i+1], CW, flag);
            }
            printf("%d\n", bs_choice.bshare);
            dump_level(nextlevel, curlevel_size<<1);
            cw.push_back(CW);
        }

        ++level;
    }

    // We don't need to store the last level

    AESkey prgkey;
    __m128i key = _mm_set_epi64x(314159265, 271828182);
    AES_128_Key_Expansion(prgkey, key);
    __m128i left, right;
    AES_ECB_encrypt(left, set_lsb(seed, 0), prgkey, aesops);
    AES_ECB_encrypt(right, set_lsb(seed, 1), prgkey, aesops);

    __m128i nleft, nright, oleft, oright;
    prg(nleft, seed, 0, aesops);
    prg(nright, seed, 1, aesops);
    prgboth(oleft, oright, seed, aesops);
    printf("left : "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&left)[15-i]); } printf("\n");
    printf("nleft: "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&nleft)[15-i]); } printf("\n");
    printf("oleft: "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&oleft)[15-i]); } printf("\n");
    printf("rght : "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&right)[15-i]); } printf("\n");
    printf("nrght: "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&nright)[15-i]); } printf("\n");
    printf("orght: "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&oright)[15-i]); } printf("\n");
}
