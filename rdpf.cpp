#include <bsd/stdlib.h> // arc4random_buf

#include "rdpf.hpp"
#include "bitutils.hpp"

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
    for(int i=0;i<16;++i) { printf("%02x", ((unsigned char *)&seed)[15-i]); } printf("\n");
    rdpf.seed = seed;
}
