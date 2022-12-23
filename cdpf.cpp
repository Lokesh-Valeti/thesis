#include <bsd/stdlib.h> // arc4random_buf
#include "bitutils.hpp"
#include "cdpf.hpp"

// Generate a pair of CDPFs with the given target value
std::tuple<CDPF,CDPF> CDPF::generate(value_t target)
{
    CDPF dpf0, dpf1;

    // Pick two random seeds
    arc4random_buf(&dpf0.seed, sizeof(dpf0.seed));
    arc4random_buf(&dpf1.seed, sizeof(dpf1.seed));
    // Ensure the flag bits (the lsb of each node) are different
    dpf0.seed = set_lsb(dpf0.seed, 0);
    dpf1.seed = set_lsb(dpf1.seed, 1);
    dpf0.whichhalf = 0;
    dpf1.whichhalf = 1;

    return std::make_tuple(dpf0, dpf1);
}

// Generate a pair of CDPFs with a random target value
std::tuple<CDPF,CDPF> CDPF::generate()
{
    value_t target;
    arc4random_buf(&target, sizeof(target));
    return generate(target);
}
