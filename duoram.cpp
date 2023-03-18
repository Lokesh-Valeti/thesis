#include "duoram.hpp"
#include "shapes.hpp"

// Assuming the memory is already sorted, do an oblivious binary
// search for the smallest index containing the value at least the
// given one.  (The answer will be the length of the Flat if all
// elements are smaller than the target.) Only available for additive
// shared databases for now.
template <>
RegAS Duoram<RegAS>::Flat::obliv_binary_search(RegAS &target)
{
    if (this->shape_size == 0) {
        RegAS zero;
        return zero;
    }
    // Create a Pad of the smallest power of 2 size strictly greater
    // than the Flat size
    address_t padsize = 1;
    nbits_t depth = 0;
    while (padsize <= this->shape_size) {
        padsize *= 2;
        ++depth;
    }
    Duoram<RegAS>::Pad P(*this, tio, yield, padsize);

    // Start in the middle
    RegAS index;
    index.set(this->tio.player() ? 0 : (1<<(depth-1))-1);
    // Invariant: index points to the last element of the left half of
    // the remaining possible range, which is of width (1<<depth).
    while (depth > 0) {
        // Obliviously read the value there
        RegAS val = P[index];
        // Compare it to the target
        CDPF cdpf = tio.cdpf(this->yield);
        auto [lt, eq, gt] = cdpf.compare(this->tio, this->yield,
            val-target, tio.aes_ops());
        if (depth > 1) {
            // If val >= target, the answer is here or to the left
            // and we should subtract 2^{depth-2} from index
            // If val < target, the answer is strictly to the right
            // and we should add 2^{depth-2} to index
            // So we unconditionally subtract 2^{depth-2} from index, and
            // add (lt)*2^{depth-1}.
            RegAS uncond;
            uncond.set(tio.player() ? 0 : address_t(1)<<(depth-2));
            RegAS cond;
            cond.set(tio.player() ? 0 : address_t(1)<<(depth-1));
            RegAS condprod;
            mpc_flagmult(this->tio, this->yield, condprod, lt, cond);
            index -= uncond;
            index += condprod;
        } else {
            // The possible range is of width 2, and we're pointing to
            // the first element of it.
            // If val >= target, the answer is here or to the left, so
            // it's here.
            // If val < target, the answer is strictly to the right
            // so add lt to index
            RegAS cond;
            cond.set(tio.player() ? 0 : 1);
            RegAS condprod;
            mpc_flagmult(this->tio, this->yield, condprod, lt, cond);
            index += condprod;
        }
        --depth;
    }

    return index;
}

