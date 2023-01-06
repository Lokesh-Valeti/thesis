#include "duoram.hpp"

// Assuming the memory is already sorted, do an oblivious binary
// search for the largest index containing the value at most the
// given one.  (The answer will be 0 if all of the memory elements
// are greate than the target.) This Flat must be a power of 2 size.
// Only available for additive shared databases for now.
template <>
RegAS Duoram<RegAS>::Flat::obliv_binary_search(RegAS &target)
{
    nbits_t depth = this->addr_size;
    // Start in the middle
    RegAS index;
    index.set(this->tio.player() ? 0 : 1<<(depth-1));
    // Invariant: index points to the first element of the right half of
    // the remaining possible range
    while (depth > 0) {
        // Obliviously read the value there
        RegAS val = operator[](index);
        // Compare it to the target
        CDPF cdpf = tio.cdpf(this->yield);
        auto [lt, eq, gt] = cdpf.compare(this->tio, this->yield,
            val-target, tio.aes_ops());
        if (depth > 1) {
            // If val > target, the answer is strictly to the left
            // and we should subtract 2^{depth-2} from index
            // If val <= target, the answer is here or to the right
            // and we should add 2^{depth-2} to index
            // So we unconditionally subtract 2^{depth-2} from index, and
            // add (lt+eq)*2^{depth-1}.
            RegAS uncond;
            uncond.set(tio.player() ? 0 : address_t(1)<<(depth-2));
            RegAS cond;
            cond.set(tio.player() ? 0 : address_t(1)<<(depth-1));
            RegAS condprod;
            RegBS le = lt ^ eq;
            mpc_flagmult(this->tio, this->yield, condprod, le, cond);
            index -= uncond;
            index += condprod;
        } else {
            // If val > target, the answer is strictly to the left
            // If val <= target, the answer is here or to the right
            // so subtract gt from index
            RegAS cond;
            cond.set(tio.player() ? 0 : 1);
            RegAS condprod;
            mpc_flagmult(this->tio, this->yield, condprod, gt, cond);
            index -= condprod;
        }
        --depth;
    }

    return index;
}

