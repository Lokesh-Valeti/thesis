#include <bsd/stdlib.h> // arc4random_buf

#include "rdpf.hpp"
#include "bitutils.hpp"
#include "aes.hpp"

// Construct a DPF of the given depth to be used for random-access
// memory reads and writes.  The DPF is construction collaboratively by
// P0 and P1, with the server P2 helping by providing various kinds of
// correlated randomness, such as MultTriples and AndTriples.
void rdpf_gen(MPCTIO &tio, yield_t &yield,
    RDPF &rdpf, nbits_t depth)
{
    int player = tio.player();

    // Choose a random seed
    DPFnode seed;
    arc4random_buf(&seed, sizeof(seed));
    // Ensure the flag bits (the lsb of each node) are different
    seed = set_lsb(seed, !!player);
    printf("seed: "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&seed)[15-i]); } printf("\n");
    rdpf.seed = seed;

    AESkey prgkey;
    __m128i key = _mm_set_epi64x(314159265, 271828182);
    AES_128_Key_Expansion(prgkey, key);
    __m128i left, right;
    AES_ECB_encrypt(left, set_lsb(seed, 0), prgkey);
    printf("left: "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&left)[15-i]); } printf("\n");
    AES_ECB_encrypt(right, set_lsb(seed, 1), prgkey);
    printf("rght: "); for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&right)[15-i]); } printf("\n");
}
