// Templated method implementations for rdpf.hpp

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

// Create a StreamEval object that will start its output at index start.
// It will wrap around to 0 when it hits 2^depth.  If use_expansion
// is true, then if the DPF has been expanded, just output values
// from that.  If use_expansion=false or if the DPF has not been
// expanded, compute the values on the fly.  If xor_offset is non-zero,
// then the outputs are actually
// DPF(start XOR xor_offset)
// DPF((start+1) XOR xor_offset)
// DPF((start+2) XOR xor_offset)
// etc.
template <typename T>
StreamEval<T>::StreamEval(const T &rdpf, address_t start,
    address_t xor_offset, size_t &aes_ops,
    bool use_expansion) : rdpf(rdpf), aes_ops(aes_ops),
    use_expansion(use_expansion), counter_xor_offset(xor_offset)
{
    depth = rdpf.depth();
    // Prevent overflow of 1<<depth
    if (depth < ADDRESS_MAX_BITS) {
        indexmask = (address_t(1)<<depth)-1;
    } else {
        indexmask = ~0;
    }
    start &= indexmask;
    counter_xor_offset &= indexmask;
    // Record that we haven't actually output the leaf for index start
    // itself yet
    nextindex = start;
    if (use_expansion && rdpf.has_expansion()) {
        // We just need to keep the counter, not compute anything
        return;
    }
    path.resize(depth);
    pathindex = start;
    path[0] = rdpf.get_seed();
    for (nbits_t i=1;i<depth;++i) {
        bool dir = !!(pathindex & (address_t(1)<<(depth-i)));
        bool xor_offset_bit =
            !!(counter_xor_offset & (address_t(1)<<(depth-i)));
        path[i] = rdpf.descend(path[i-1], i-1,
            dir ^ xor_offset_bit, aes_ops);
    }
}

template <typename T>
typename T::LeafNode StreamEval<T>::next()
{
    if (use_expansion && rdpf.has_expansion()) {
        // Just use the precomputed values
        typename T::LeafNode leaf =
            rdpf.get_expansion(nextindex ^ counter_xor_offset);
        nextindex = (nextindex + 1) & indexmask;
        return leaf;
    }
    // Invariant: in the first call to next(), nextindex = pathindex.
    // Otherwise, nextindex = pathindex+1.
    // Get the XOR of nextindex and pathindex, and strip the low bit.
    // If nextindex and pathindex are equal, or pathindex is even
    // and nextindex is the consecutive odd number, index_xor will be 0,
    // indicating that we don't have to update the path, but just
    // compute the appropriate leaf given by the low bit of nextindex.
    //
    // Otherwise, say for example pathindex is 010010111 and nextindex
    // is 010011000.  Then their XOR is 000001111, and stripping the low
    // bit yields 000001110, so how_many_1_bits will be 3.
    // That indicates (typically) that path[depth-3] was a left child,
    // and now we need to change it to a right child by descending right
    // from path[depth-4], and then filling the path after that with
    // left children.
    //
    // When we wrap around, however, index_xor will be 111111110 (after
    // we strip the low bit), and how_many_1_bits will be depth-1, but
    // the new top child (of the root seed) we have to compute will be a
    // left, not a right, child.
    uint64_t index_xor = (nextindex ^ pathindex) & ~1;
    nbits_t how_many_1_bits = __builtin_popcount(index_xor);
    if (how_many_1_bits > 0) {
        // This will almost always be 1, unless we've just wrapped
        // around from the right subtree back to the left, in which case
        // it will be 0.
        bool top_changed_bit =
            !!(nextindex & (address_t(1) << how_many_1_bits));
        bool xor_offset_bit =
            !!(counter_xor_offset & (address_t(1) << how_many_1_bits));
        path[depth-how_many_1_bits] =
            rdpf.descend(path[depth-how_many_1_bits-1],
                depth-how_many_1_bits-1,
                top_changed_bit ^ xor_offset_bit, aes_ops);
        for (nbits_t i = depth-how_many_1_bits; i < depth-1; ++i) {
            bool xor_offset_bit =
                !!(counter_xor_offset & (address_t(1) << (depth-i-1)));
            path[i+1] = rdpf.descend(path[i], i, xor_offset_bit, aes_ops);
        }
    }
    bool xor_offset_bit = counter_xor_offset & 1;
    typename T::LeafNode leaf = rdpf.descend_to_leaf(path[depth-1], depth-1,
        (nextindex & 1) ^ xor_offset_bit, aes_ops);
    pathindex = nextindex;
    nextindex = (nextindex + 1) & indexmask;
    return leaf;
}

// Run the parallel evaluator.  The type V is the type of the
// accumulator; init should be the "zero" value of the accumulator.
// The type W (process) is a lambda type with the signature
// (int, address_t, const T::node &) -> V
// which will be called like this for each i from 0 to num_evals-1,
// across num_thread threads:
// value_i = process(t, i, DPF((start+i) XOR xor_offset))
// t is the thread number (0 <= t < num_threads).
// The resulting num_evals values will be combined using V's +=
// operator, first accumulating the values within each thread
// (starting with the init value), and then accumulating the totals
// from each thread together (again starting with the init value):
//
// total = init
// for each thread t:
//     accum_t = init
//     for each accum_i generated by thread t:
//         accum_t += value_i
//     total += accum_t
template <typename T> template <typename V, typename W>
inline V ParallelEval<T>::reduce(V init, W process)
{
    size_t thread_aes_ops[num_threads];
    V accums[num_threads];
    boost::asio::thread_pool pool(num_threads);
    address_t threadstart = start;
    address_t threadchunk = num_evals / num_threads;
    address_t threadextra = num_evals % num_threads;
    nbits_t depth = rdpf.depth();
    address_t indexmask = (depth < ADDRESS_MAX_BITS ?
        ((address_t(1)<<depth)-1) : ~0);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        address_t threadsize = threadchunk + (address_t(thread_num) < threadextra);
        boost::asio::post(pool,
            [this, &init, &thread_aes_ops, &accums, &process,
                    thread_num, threadstart, threadsize, indexmask] {
                size_t local_aes_ops = 0;
                auto ev = StreamEval(rdpf, (start+threadstart)&indexmask,
                    xor_offset, local_aes_ops);
                V accum = init;
                for (address_t x=0;x<threadsize;++x) {
                    typename T::LeafNode leaf = ev.next();
                    accum += process(thread_num,
                        (threadstart+x)&indexmask, leaf);
                }
                accums[thread_num] = accum;
                thread_aes_ops[thread_num] = local_aes_ops;
            });
        threadstart = (threadstart + threadsize) & indexmask;
    }
    pool.join();
    V total = init;
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        total += accums[thread_num];
        aes_ops += thread_aes_ops[thread_num];
    }
    return total;
}

// Descend from a node at depth parentdepth to one of its leaf children
// whichchild = 0: left child
// whichchild = 1: right child
//
// Cost: 1 AES operation
template <nbits_t WIDTH>
inline typename RDPF<WIDTH>::LeafNode RDPF<WIDTH>::descend_to_leaf(
    const DPFnode &parent, nbits_t parentdepth, bit_t whichchild,
    size_t &aes_ops) const
{
    typename RDPF<WIDTH>::LeafNode prgout;
    bool flag = get_lsb(parent);
    prgleaf(prgout, parent, whichchild, aes_ops);
    if (flag) {
        LeafNode CW = li[0].leaf_cw;
        LeafNode CWR = CW;
        bit_t cfbit = !!(leaf_cfbits &
            (value_t(1)<<(maxdepth-parentdepth)));
        CWR[0] ^= lsb128_mask[cfbit];
        prgout ^= (whichchild ? CWR : CW);
    }
    return prgout;
}

// I/O for RDPFs

template <typename T, nbits_t WIDTH>
T& operator>>(T &is, RDPF<WIDTH> &rdpf)
{
    is.read((char *)&rdpf.seed, sizeof(rdpf.seed));
    rdpf.whichhalf = get_lsb(rdpf.seed);
    uint8_t depth;
    // Add 64 to depth to indicate an expanded RDPF
    is.read((char *)&depth, sizeof(depth));
    bool read_expanded = false;
    if (depth > 64) {
        read_expanded = true;
        depth -= 64;
    }
    assert(depth <= ADDRESS_MAX_BITS);
    rdpf.cw.clear();
    for (uint8_t i=0; i<depth; ++i) {
        DPFnode cw;
        is.read((char *)&cw, sizeof(cw));
        rdpf.cw.push_back(cw);
    }
    if (read_expanded) {
        rdpf.expansion.resize(1<<depth);
        is.read((char *)rdpf.expansion.data(),
            sizeof(rdpf.expansion[0])<<depth);
    }
    value_t cfbits = 0;
    is.read((char *)&cfbits, BITBYTES(depth));
    rdpf.cfbits = cfbits;
    rdpf.li.resize(1);
    is.read((char *)&rdpf.li[0].unit_sum_inverse,
        sizeof(rdpf.li[0].unit_sum_inverse));
    is.read((char *)&rdpf.li[0].scaled_sum,
        sizeof(rdpf.li[0].scaled_sum));
    is.read((char *)&rdpf.li[0].scaled_xor,
        sizeof(rdpf.li[0].scaled_xor));

    return is;
}

// Write the DPF to the output stream.  If expanded=true, then include
// the expansion _if_ the DPF is itself already expanded.  You can use
// this to write DPFs to files.
template <typename T, nbits_t WIDTH>
T& write_maybe_expanded(T &os, const RDPF<WIDTH> &rdpf,
    bool expanded = true)
{
    os.write((const char *)&rdpf.seed, sizeof(rdpf.seed));
    uint8_t depth = rdpf.cw.size();
    assert(depth <= ADDRESS_MAX_BITS);
    // If we're writing an expansion, add 64 to depth
    uint8_t expanded_depth = depth;
    bool write_expansion = false;
    if (expanded && rdpf.expansion.size() == (size_t(1)<<depth)) {
        write_expansion = true;
        expanded_depth += 64;
    }
    os.write((const char *)&expanded_depth, sizeof(expanded_depth));
    for (uint8_t i=0; i<depth; ++i) {
        os.write((const char *)&rdpf.cw[i], sizeof(rdpf.cw[i]));
    }
    if (write_expansion) {
        os.write((const char *)rdpf.expansion.data(),
            sizeof(rdpf.expansion[0])<<depth);
    }
    os.write((const char *)&rdpf.cfbits, BITBYTES(depth));
    os.write((const char *)&rdpf.li[0].unit_sum_inverse,
        sizeof(rdpf.li[0].unit_sum_inverse));
    os.write((const char *)&rdpf.li[0].scaled_sum,
        sizeof(rdpf.li[0].scaled_sum));
    os.write((const char *)&rdpf.li[0].scaled_xor,
        sizeof(rdpf.li[0].scaled_xor));

    return os;
}

// The ordinary << version never writes the expansion, since this is
// what we use to send DPFs over the network.
template <typename T, nbits_t WIDTH>
T& operator<<(T &os, const RDPF<WIDTH> &rdpf)
{
    return write_maybe_expanded(os, rdpf, false);
}

// I/O for RDPF Triples

// We never write RDPFTriples over the network, so always write
// the DPF expansions if they're available.
template <typename T, nbits_t WIDTH>
T& operator<<(T &os, const RDPFTriple<WIDTH> &rdpftrip)
{
    write_maybe_expanded(os, rdpftrip.dpf[0], true);
    write_maybe_expanded(os, rdpftrip.dpf[1], true);
    write_maybe_expanded(os, rdpftrip.dpf[2], true);
    nbits_t depth = rdpftrip.dpf[0].depth();
    os.write((const char *)&rdpftrip.as_target.ashare, BITBYTES(depth));
    os.write((const char *)&rdpftrip.xs_target.xshare, BITBYTES(depth));
    return os;
}

template <typename T, nbits_t WIDTH>
T& operator>>(T &is, RDPFTriple<WIDTH> &rdpftrip)
{
    is >> rdpftrip.dpf[0] >> rdpftrip.dpf[1] >> rdpftrip.dpf[2];
    nbits_t depth = rdpftrip.dpf[0].depth();
    rdpftrip.as_target.ashare = 0;
    is.read((char *)&rdpftrip.as_target.ashare, BITBYTES(depth));
    rdpftrip.xs_target.xshare = 0;
    is.read((char *)&rdpftrip.xs_target.xshare, BITBYTES(depth));
    return is;
}

// I/O for RDPF Pairs

// We never write RDPFPairs over the network, so always write
// the DPF expansions if they're available.
template <typename T, nbits_t WIDTH>
T& operator<<(T &os, const RDPFPair<WIDTH> &rdpfpair)
{
    write_maybe_expanded(os, rdpfpair.dpf[0], true);
    write_maybe_expanded(os, rdpfpair.dpf[1], true);
    return os;
}

template <typename T, nbits_t WIDTH>
T& operator>>(T &is, RDPFPair<WIDTH> &rdpfpair)
{
    is >> rdpfpair.dpf[0] >> rdpfpair.dpf[1];
    return is;
}

// Construct a DPF with the given (XOR-shared) target location, and
// of the given depth, to be used for random-access memory reads and
// writes.  The DPF is construction collaboratively by P0 and P1,
// with the server P2 helping by providing various kinds of
// correlated randomness, such as MultTriples and AndTriples.
//
// This algorithm is based on Appendix C from the Duoram paper, with a
// small optimization noted below.
template <nbits_t WIDTH>
RDPF<WIDTH>::RDPF(MPCTIO &tio, yield_t &yield,
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
    maxdepth = depth;
    curdepth = depth;

    // The root level is just the seed
    nbits_t level = 0;
    DPFnode *curlevel = NULL;
    DPFnode *nextlevel = new DPFnode[1];
    nextlevel[0] = seed;

    li.resize(1);

    // Construct each intermediate level
    while(level < depth) {
        if (player < 2) {
            delete[] curlevel;
            curlevel = nextlevel;
            if (save_expansion && level == depth-1) {
                expansion.resize(1<<depth);
                nextlevel = (DPFnode *)expansion.data();
            } else {
                nextlevel = new DPFnode[1<<(level+1)];
            }
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
#ifdef RDPF_MTGEN_TIMING_1
            if (player == 0) {
                mtgen_timetest_1(level, 0, (1<<23)>>level, curlevel,
                    nextlevel, aes_ops);
                size_t niters = 2048;
                if (level > 8) niters = (1<<20)>>level;
                for(int t=1;t<=8;++t) {
                    mtgen_timetest_1(level, t, niters, curlevel,
                        nextlevel, aes_ops);
                }
                mtgen_timetest_1(level, 0, (1<<23)>>level, curlevel,
                    nextlevel, aes_ops);
            }
#endif
            // Using the timing results gathered above, decide whether
            // to multithread, and if so, how many threads to use.
            // tio.cpu_nthreads() is the maximum number we have
            // available.
            int max_nthreads = tio.cpu_nthreads();
            if (max_nthreads == 1 || level < 19) {
                // No threading
                size_t laes_ops = 0;
                for(size_t i=0;i<curlevel_size;++i) {
                    DPFnode lchild, rchild;
                    prgboth(lchild, rchild, curlevel[i], laes_ops);
                    L = (L ^ lchild);
                    R = (R ^ rchild);
                    nextlevel[2*i] = lchild;
                    nextlevel[2*i+1] = rchild;
                }
                aes_ops += laes_ops;
            } else {
                size_t curlevel_size = size_t(1)<<level;
                int nthreads =
                    int(ceil(sqrt(double(curlevel_size/6000))));
                if (nthreads > max_nthreads) {
                    nthreads = max_nthreads;
                }
                DPFnode tL[nthreads];
                DPFnode tR[nthreads];
                size_t taes_ops[nthreads];
                size_t threadstart = 0;
                size_t threadchunk = curlevel_size / nthreads;
                size_t threadextra = curlevel_size % nthreads;
                boost::asio::thread_pool pool(nthreads);
                for (int t=0;t<nthreads;++t) {
                    size_t threadsize = threadchunk + (size_t(t) < threadextra);
                    size_t threadend = threadstart + threadsize;
                    boost::asio::post(pool,
                        [t, &tL, &tR, &taes_ops, threadstart, threadend,
                        &curlevel, &nextlevel] {
                            DPFnode L = _mm_setzero_si128();
                            DPFnode R = _mm_setzero_si128();
                            size_t aes_ops = 0;
                            for(size_t i=threadstart;i<threadend;++i) {
                                DPFnode lchild, rchild;
                                prgboth(lchild, rchild, curlevel[i], aes_ops);
                                L = (L ^ lchild);
                                R = (R ^ rchild);
                                nextlevel[2*i] = lchild;
                                nextlevel[2*i+1] = rchild;
                            }
                            tL[t] = L;
                            tR[t] = R;
                            taes_ops[t] = aes_ops;
                        });
                    threadstart = threadend;
                }
                pool.join();
                for (int t=0;t<nthreads;++t) {
                    L ^= tL[t];
                    R ^= tR[t];
                    aes_ops += taes_ops[t];
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
        run_coroutines(yield,
            [this, &tio, &our_parity_bit, &peer_parity_bit](yield_t &yield) {
                tio.queue_peer(&our_parity_bit, 1);
                yield();
                uint8_t peer_parity_byte;
                tio.recv_peer(&peer_parity_byte, 1);
                peer_parity_bit = peer_parity_byte & 1;
            },
            [this, &tio, &CW, &L, &R, &bs_choice, &our_parity](yield_t &yield) {
                mpc_reconstruct_choice(tio, yield, CW, bs_choice,
                    (R ^ our_parity), L);
            });
        bool parity_bit = our_parity_bit ^ peer_parity_bit;
        cfbits |= (value_t(parity_bit)<<level);
        DPFnode CWR = CW ^ lsb128_mask[parity_bit];
        if (player < 2) {
            // The timing of each iteration of the inner loop is
            // comparable to the above, so just use the same
            // computations.  All of this could be tuned, of course.

            if (level < depth-1) {
                // Using the timing results gathered above, decide whether
                // to multithread, and if so, how many threads to use.
                // tio.cpu_nthreads() is the maximum number we have
                // available.
                int max_nthreads = tio.cpu_nthreads();
                if (max_nthreads == 1 || level < 19) {
                    // No threading
                    for(size_t i=0;i<curlevel_size;++i) {
                        bool flag = get_lsb(curlevel[i]);
                        nextlevel[2*i] = xor_if(nextlevel[2*i], CW, flag);
                        nextlevel[2*i+1] = xor_if(nextlevel[2*i+1], CWR, flag);
                    }
                } else {
                    int nthreads =
                        int(ceil(sqrt(double(curlevel_size/6000))));
                    if (nthreads > max_nthreads) {
                        nthreads = max_nthreads;
                    }
                    size_t threadstart = 0;
                    size_t threadchunk = curlevel_size / nthreads;
                    size_t threadextra = curlevel_size % nthreads;
                    boost::asio::thread_pool pool(nthreads);
                    for (int t=0;t<nthreads;++t) {
                        size_t threadsize = threadchunk + (size_t(t) < threadextra);
                        size_t threadend = threadstart + threadsize;
                        boost::asio::post(pool, [CW, CWR, threadstart, threadend,
                            &curlevel, &nextlevel] {
                                for(size_t i=threadstart;i<threadend;++i) {
                                    bool flag = get_lsb(curlevel[i]);
                                    nextlevel[2*i] = xor_if(nextlevel[2*i], CW, flag);
                                    nextlevel[2*i+1] = xor_if(nextlevel[2*i+1], CWR, flag);
                                }
                        });
                        threadstart = threadend;
                    }
                    pool.join();
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
                // Using the timing results gathered above, decide whether
                // to multithread, and if so, how many threads to use.
                // tio.cpu_nthreads() is the maximum number we have
                // available.
                int max_nthreads = tio.cpu_nthreads();
                if (max_nthreads == 1 || level < 19) {
                    // No threading
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
                } else {
                    int nthreads =
                        int(ceil(sqrt(double(curlevel_size/6000))));
                    if (nthreads > max_nthreads) {
                        nthreads = max_nthreads;
                    }
                    value_t tlow_sum[nthreads];
                    value_t thigh_sum[nthreads];
                    value_t thigh_xor[nthreads];
                    size_t threadstart = 0;
                    size_t threadchunk = curlevel_size / nthreads;
                    size_t threadextra = curlevel_size % nthreads;
                    boost::asio::thread_pool pool(nthreads);
                    for (int t=0;t<nthreads;++t) {
                        size_t threadsize = threadchunk + (size_t(t) < threadextra);
                        size_t threadend = threadstart + threadsize;
                        boost::asio::post(pool,
                            [t, &tlow_sum, &thigh_sum, &thigh_xor, threadstart, threadend,
                            &curlevel, &nextlevel, CW, CWR, save_expansion] {
                                value_t low_sum = 0;
                                value_t high_sum = 0;
                                value_t high_xor = 0;
                                for(size_t i=threadstart;i<threadend;++i) {
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
                                tlow_sum[t] = low_sum;
                                thigh_sum[t] = high_sum;
                                thigh_xor[t] = high_xor;
                            });
                        threadstart = threadend;
                    }
                    pool.join();
                    for (int t=0;t<nthreads;++t) {
                        low_sum += tlow_sum[t];
                        high_sum += thigh_sum[t];
                        high_xor ^= thigh_xor[t];
                    }
                }
                if (player == 1) {
                    low_sum = -low_sum;
                    high_sum = -high_sum;
                }
                li[0].scaled_sum[0].ashare = high_sum;
                li[0].scaled_xor[0].xshare = high_xor;
                // Exchange low_sum and add them up
                tio.queue_peer(&low_sum, sizeof(low_sum));
                yield();
                value_t peer_low_sum;
                tio.recv_peer(&peer_low_sum, sizeof(peer_low_sum));
                low_sum += peer_low_sum;
                // The low_sum had better be odd
                assert(low_sum & 1);
                li[0].unit_sum_inverse = inverse_value_t(low_sum);
            }
            cw.push_back(CW);
        } else if (level == depth-1) {
            yield();
        }

        ++level;
    }

    delete[] curlevel;
    if (!save_expansion || player == 2) {
        delete[] nextlevel;
    }
}

// Get the leaf node for the given input
template <nbits_t WIDTH>
typename RDPF<WIDTH>::LeafNode
    RDPF<WIDTH>::leaf(address_t input, size_t &aes_ops) const
{
    // If we have a precomputed expansion, just use it
    if (expansion.size()) {
        return expansion[input];
    }

    DPFnode node = seed;
    for (nbits_t d=0;d<curdepth-1;++d) {
        bit_t dir = !!(input & (address_t(1)<<(curdepth-d-1)));
        node = descend(node, d, dir, aes_ops);
    }
    bit_t dir = (input & 1);
    return descend_to_leaf(node, curdepth, dir, aes_ops);
}

// Expand the DPF if it's not already expanded
//
// This routine is slightly more efficient than repeatedly calling
// StreamEval::next(), but it uses a lot more memory.
template <nbits_t WIDTH>
void RDPF<WIDTH>::expand(size_t &aes_ops)
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
        path[i] = descend(path[i-1], i-1, 0, aes_ops);
    }
    expansion[index++] = descend_to_leaf(path[depth-1], depth-1, 0, aes_ops);
    expansion[index++] = descend_to_leaf(path[depth-1], depth-1, 1, aes_ops);
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
                depth-how_many_1_bits-1, 1, aes_ops);
        for (nbits_t i = depth-how_many_1_bits; i < depth-1; ++i) {
            path[i+1] = descend(path[i], i, 0, aes_ops);
        }
        lastindex = index;
        expansion[index++] = descend_to_leaf(path[depth-1], depth-1, 0, aes_ops);
        expansion[index++] = descend_to_leaf(path[depth-1], depth-1, 1, aes_ops);
    }

    delete[] path;
}

// Construct three RDPFs of the given depth all with the same randomly
// generated target index.
template <nbits_t WIDTH>
RDPFTriple<WIDTH>::RDPFTriple(MPCTIO &tio, yield_t &yield,
    nbits_t depth, bool save_expansion)
{
    // Pick a random XOR share of the target
    xs_target.randomize(depth);

    // Now create three RDPFs with that target, and also convert the XOR
    // shares of the target to additive shares
    std::vector<coro_t> coroutines;
    for (int i=0;i<3;++i) {
        coroutines.emplace_back(
            [this, &tio, depth, i, save_expansion](yield_t &yield) {
                dpf[i] = RDPF<WIDTH>(tio, yield, xs_target, depth,
                    save_expansion);
            });
    }
    coroutines.emplace_back(
        [this, &tio, depth](yield_t &yield) {
            mpc_xs_to_as(tio, yield, as_target, xs_target, depth, false);
        });
    run_coroutines(yield, coroutines);
}

template <nbits_t WIDTH>
typename RDPFTriple<WIDTH>::node RDPFTriple<WIDTH>::descend(
    const RDPFTriple<WIDTH>::node &parent,
    nbits_t parentdepth, bit_t whichchild,
    size_t &aes_ops) const
{
    auto [P0, P1, P2] = parent;
    DPFnode C0, C1, C2;
    C0 = dpf[0].descend(P0, parentdepth, whichchild, aes_ops);
    C1 = dpf[1].descend(P1, parentdepth, whichchild, aes_ops);
    C2 = dpf[2].descend(P2, parentdepth, whichchild, aes_ops);
    return std::make_tuple(C0,C1,C2);
}

template <nbits_t WIDTH>
typename RDPFTriple<WIDTH>::LeafNode RDPFTriple<WIDTH>::descend_to_leaf(
    const RDPFTriple<WIDTH>::node &parent,
    nbits_t parentdepth, bit_t whichchild,
    size_t &aes_ops) const
{
    auto [P0, P1, P2] = parent;
    typename RDPF<WIDTH>::LeafNode C0, C1, C2;
    C0 = dpf[0].descend_to_leaf(P0, parentdepth, whichchild, aes_ops);
    C1 = dpf[1].descend_to_leaf(P1, parentdepth, whichchild, aes_ops);
    C2 = dpf[2].descend_to_leaf(P2, parentdepth, whichchild, aes_ops);
    return std::make_tuple(C0,C1,C2);
}

template <nbits_t WIDTH>
typename RDPFPair<WIDTH>::node RDPFPair<WIDTH>::descend(
    const RDPFPair<WIDTH>::node &parent,
    nbits_t parentdepth, bit_t whichchild,
    size_t &aes_ops) const
{
    auto [P0, P1] = parent;
    DPFnode C0, C1;
    C0 = dpf[0].descend(P0, parentdepth, whichchild, aes_ops);
    C1 = dpf[1].descend(P1, parentdepth, whichchild, aes_ops);
    return std::make_tuple(C0,C1);
}

template <nbits_t WIDTH>
typename RDPFPair<WIDTH>::LeafNode RDPFPair<WIDTH>::descend_to_leaf(
    const RDPFPair<WIDTH>::node &parent,
    nbits_t parentdepth, bit_t whichchild,
    size_t &aes_ops) const
{
    auto [P0, P1] = parent;
    typename RDPF<WIDTH>::LeafNode C0, C1;
    C0 = dpf[0].descend_to_leaf(P0, parentdepth, whichchild, aes_ops);
    C1 = dpf[1].descend_to_leaf(P1, parentdepth, whichchild, aes_ops);
    return std::make_tuple(C0,C1);
}
