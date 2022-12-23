#include <bsd/stdlib.h> // arc4random_buf
#include "bitutils.hpp"
#include "cdpf.hpp"

// Generate a pair of CDPFs with the given target value
std::tuple<CDPF,CDPF> CDPF::generate(value_t target, size_t &aes_ops)
{
    CDPF dpf0, dpf1;
    nbits_t depth = VALUE_BITS - 7;

    // Pick two random seeds
    arc4random_buf(&dpf0.seed, sizeof(dpf0.seed));
    arc4random_buf(&dpf1.seed, sizeof(dpf1.seed));
    // Ensure the flag bits (the lsb of each node) are different
    dpf0.seed = set_lsb(dpf0.seed, 0);
    dpf1.seed = set_lsb(dpf1.seed, 1);
    dpf0.whichhalf = 0;
    dpf1.whichhalf = 1;
    dpf0.cfbits = 0;
    dpf1.cfbits = 0;
    dpf0.as_target.randomize();
    dpf1.as_target.ashare = target - dpf0.as_target.ashare;
    dpf0.xs_target.randomize();
    dpf1.xs_target.xshare = target ^ dpf0.xs_target.xshare;

    // The current node in each CDPF as we descend the tree.  The
    // invariant is that cur0 and cur1 are the nodes on the path to the
    // target at level curlevel.  They will necessarily be different,
    // and indeed must have different flag (low) bits.
    DPFnode cur0 = dpf0.seed;
    DPFnode cur1 = dpf1.seed;
    nbits_t curlevel = 0;

    while(curlevel < depth) {
        // Construct the two (uncorrected) children of each node
        DPFnode left0, right0, left1, right1;
        prgboth(left0, right0, cur0, aes_ops);
        prgboth(left1, right1, cur1, aes_ops);

        // Which way lies the target?
        bool targetdir = !!(target & (value_t(1)<<(depth-curlevel-1)));
        DPFnode CW;
        bool cfbit = !get_lsb(left0 ^ left1 ^ right0 ^ right1);
        bool flag0 = get_lsb(cur0);
        bool flag1 = get_lsb(cur1);
        // The last level is special
        if (curlevel < depth-1) {
            if (targetdir == 0) {
                // The target is to the left, so make the correction word
                // and bit make the right children the same and the left
                // children have different flag bits.

                // Recall that descend will apply (only for the party whose
                // current node (cur0 or cur1) has the flag bit set, for
                // which exactly one of the two will) CW to both children,
                // and cfbit to the flag bit of the right child.
                CW = right0 ^ right1 ^ lsb128_mask[cfbit];

                // Compute the current nodes for the next level
                // Exactly one of these two XORs will fire, so afterwards,
                // cur0 ^ cur1 = left0 ^ left1 ^ CW, which will have low bit
                // 1 by the definition of cfbit.
                cur0 = xor_if(left0, CW, flag0);
                cur1 = xor_if(left1, CW, flag1);
            } else {
                // The target is to the right, so make the correction word
                // and bit make the left children the same and the right
                // children have different flag bits.
                CW = left0 ^ left1;

                // Compute the current nodes for the next level
                // Exactly one of these two XORs will fire, so similar to
                // the above, afterwards, cur0 ^ cur1 = right0 ^ right1 ^ CWR,
                // which will have low bit 1.
                DPFnode CWR = CW ^ lsb128_mask[cfbit];
                cur0 = xor_if(right0, CWR, flag0);
                cur1 = xor_if(right1, CWR, flag1);
            }
        } else {
            // We're at the last level before the leaves.  We still want
            // the children not in the direction of targetdir to end up
            // the same, but now we want the child in the direction of
            // targetdir to also end up the same, except for the single
            // target bit.  Importantly, the low bit (the flag bit in
            // all other nodes) is not special, and will in fact usually
            // end up the same for the two DPFs (unless the target bit
            // happens to be the low bit of the word; i.e., the low 7
            // bits of target are all 0).

            // This will be a 128-bit word with a single bit set, in
            // position (target & 0x7f).
            uint8_t loc = (target & 0x7f);
            DPFnode target_set_bit = _mm_set_epi64x(
                loc >= 64 ? (uint64_t(1)<<(loc-64)) : 0,
                loc >= 64 ? 0 : (uint64_t(1)<<loc));

            if (targetdir == 0) {
                // We want the right children to be the same, and the
                // left children to be the same except for the target
                // bit.
                // Remember for exactly one of the two parties, CW will
                // be applied to the left and CWR will be applied to the
                // right.
                CW = left0 ^ left1 ^ target_set_bit;
                DPFnode CWR = right0 ^ right1;
                dpf0.leaf_cwr = CWR;
                dpf1.leaf_cwr = CWR;
            } else {
                // We want the left children to be the same, and the
                // right children to be the same except for the target
                // bit.
                // Remember for exactly one of the two parties, CW will
                // be applied to the left and CWR will be applied to the
                // right.
                CW = left0 ^ left1;
                DPFnode CWR = right0 ^ right1 ^ target_set_bit;
                dpf0.leaf_cwr = CWR;
                dpf1.leaf_cwr = CWR;
            }
        }
        dpf0.cw.push_back(CW);
        dpf1.cw.push_back(CW);
        dpf0.cfbits |= (value_t(cfbit)<<curlevel);
        dpf1.cfbits |= (value_t(cfbit)<<curlevel);
        ++curlevel;
    }

    return std::make_tuple(dpf0, dpf1);
}

// Generate a pair of CDPFs with a random target value
std::tuple<CDPF,CDPF> CDPF::generate(size_t &aes_ops)
{
    value_t target;
    arc4random_buf(&target, sizeof(target));
    return generate(target, aes_ops);
}
