// Templated method implementations for rdpf.hpp

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
typename T::node StreamEval<T>::next()
{
    if (use_expansion && rdpf.has_expansion()) {
        // Just use the precomputed values
        typename T::node leaf =
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
    typename T::node leaf = rdpf.descend(path[depth-1], depth-1,
        (nextindex & 1) ^ xor_offset_bit, aes_ops);
    pathindex = nextindex;
    nextindex = (nextindex + 1) & indexmask;
    return leaf;
}

// Run the parallel evaluator.  The type V is the type of the
// accumulator; init should be the "zero" value of the accumulator.
// The type W (process) is a lambda type with the signature
// (const ParallelEval &, int, address_t, T::node) -> V
// which will be called like this for each i from 0 to num_evals-1,
// across num_thread threads:
// value_i = process(*this, t, i, DPF((start+i) XOR xor_offset))
// t is the thread number (0 <= t < num_threads).
// The type X (accumulate) is a lambda type with the signature
// (const ParallelEval &, V &, const V &)
// which will be called to combine the num_evals values of accum,
// first accumulating the values within each thread (starting with
// the init value), and then accumulating the totals from each
// thread together (again starting with the init value):
//
// total = init
// for each thread t:
//     accum_t = init
//     for each accum_i generated by thread t:
//         accumulate(*this, acccum_t, value_i)
//     accumulate(*this, total, accum_t)
template <typename T> template <typename V, typename W, typename X>
inline V ParallelEval<T>::reduce(V init, W process, X accumulate)
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
                    &accumulate, thread_num, threadstart, threadsize,
                    indexmask] {
                size_t local_aes_ops = 0;
                auto ev = StreamEval(rdpf, (start+threadstart)&indexmask,
                    xor_offset, local_aes_ops);
                V accum = init;
                for (address_t x=0;x<threadsize;++x) {
                    typename T::node leaf = ev.next();
                    V value = process(*this, thread_num,
                        (threadstart+x)&indexmask, leaf);
                    accumulate(*this, accum, value);
                }
                accums[thread_num] = accum;
                thread_aes_ops[thread_num] = local_aes_ops;
            });
        threadstart = (threadstart + threadsize) & indexmask;
    }
    pool.join();
    V total = init;
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        accumulate(*this, total, accums[thread_num]);
        aes_ops += thread_aes_ops[thread_num];
    }
    return total;
}

// Additive share of the target index
template <>
inline RegAS RDPFTriple::target<RegAS>() const {
    return as_target;
}

// XOR share of the target index
template <>
inline RegXS RDPFTriple::target<RegXS>() const {
    return xs_target;
}

// Additive share of the scaling value M_as such that the high words
// of the leaf values for P0 and P1 add to M_as * e_{target}
template <>
inline std::tuple<RegAS,RegAS,RegAS> RDPFTriple::scaled_value<RegAS>() const {
    return std::make_tuple(dpf[0].scaled_sum, dpf[1].scaled_sum,
        dpf[2].scaled_sum);
}

// XOR share of the scaling value M_xs such that the high words
// of the leaf values for P0 and P1 XOR to M_xs * e_{target}
template <>
inline std::tuple<RegXS,RegXS,RegXS> RDPFTriple::scaled_value<RegXS>() const {
    return std::make_tuple(dpf[0].scaled_xor, dpf[1].scaled_xor,
        dpf[2].scaled_xor);
}

// Get the bit-shared unit vector entry from the leaf node
template <>
inline std::tuple<RegXS,RegXS,RegXS> RDPFTriple::unit<RegXS>(node leaf) const {
    return std::make_tuple(
        dpf[0].unit_bs(std::get<0>(leaf)),
        dpf[1].unit_bs(std::get<1>(leaf)),
        dpf[2].unit_bs(std::get<2>(leaf)));
}

// Get the additive-shared unit vector entry from the leaf node
template <>
inline std::tuple<RegAS,RegAS,RegAS> RDPFTriple::unit<RegAS>(node leaf) const {
    return std::make_tuple(
        dpf[0].unit_as(std::get<0>(leaf)),
        dpf[1].unit_as(std::get<1>(leaf)),
        dpf[2].unit_as(std::get<2>(leaf)));
}

// Get the XOR-shared scaled vector entry from the leaf ndoe
template <>
inline std::tuple<RegXS,RegXS,RegXS> RDPFTriple::scaled<RegXS>(node leaf) const {
    return std::make_tuple(
        dpf[0].scaled_xs(std::get<0>(leaf)),
        dpf[1].scaled_xs(std::get<1>(leaf)),
        dpf[2].scaled_xs(std::get<2>(leaf)));
}

// Get the additive-shared scaled vector entry from the leaf node
template <>
inline std::tuple<RegAS,RegAS,RegAS> RDPFTriple::scaled<RegAS>(node leaf) const {
    return std::make_tuple(
        dpf[0].scaled_as(std::get<0>(leaf)),
        dpf[1].scaled_as(std::get<1>(leaf)),
        dpf[2].scaled_as(std::get<2>(leaf)));
}

// Additive share of the scaling value M_as such that the high words
// of the leaf values for P0 and P1 add to M_as * e_{target}
template <>
inline std::tuple<RegAS,RegAS> RDPFPair::scaled_value<RegAS>() const {
    return std::make_tuple(dpf[0].scaled_sum, dpf[1].scaled_sum);
}

// XOR share of the scaling value M_xs such that the high words
// of the leaf values for P0 and P1 XOR to M_xs * e_{target}
template <>
inline std::tuple<RegXS,RegXS> RDPFPair::scaled_value<RegXS>() const {
    return std::make_tuple(dpf[0].scaled_xor, dpf[1].scaled_xor);
}

// Get the bit-shared unit vector entry from the leaf node
template <>
inline std::tuple<RegXS,RegXS> RDPFPair::unit<RegXS>(node leaf) const {
    return std::make_tuple(
        dpf[0].unit_bs(std::get<0>(leaf)),
        dpf[1].unit_bs(std::get<1>(leaf)));
}

// Get the additive-shared unit vector entry from the leaf node
template <>
inline std::tuple<RegAS,RegAS> RDPFPair::unit<RegAS>(node leaf) const {
    return std::make_tuple(
        dpf[0].unit_as(std::get<0>(leaf)),
        dpf[1].unit_as(std::get<1>(leaf)));
}

// Get the XOR-shared scaled vector entry from the leaf ndoe
template <>
inline std::tuple<RegXS,RegXS> RDPFPair::scaled<RegXS>(node leaf) const {
    return std::make_tuple(
        dpf[0].scaled_xs(std::get<0>(leaf)),
        dpf[1].scaled_xs(std::get<1>(leaf)));
}

// Get the additive-shared scaled vector entry from the leaf node
template <>
inline std::tuple<RegAS,RegAS> RDPFPair::scaled<RegAS>(node leaf) const {
    return std::make_tuple(
        dpf[0].scaled_as(std::get<0>(leaf)),
        dpf[1].scaled_as(std::get<1>(leaf)));
}

// I/O for RDPFs

template <typename T>
T& operator>>(T &is, RDPF &rdpf)
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
    is.read((char *)&rdpf.unit_sum_inverse, sizeof(rdpf.unit_sum_inverse));
    is.read((char *)&rdpf.scaled_sum, sizeof(rdpf.scaled_sum));
    is.read((char *)&rdpf.scaled_xor, sizeof(rdpf.scaled_xor));

    return is;
}

// Write the DPF to the output stream.  If expanded=true, then include
// the expansion _if_ the DPF is itself already expanded.  You can use
// this to write DPFs to files.
template <typename T>
T& write_maybe_expanded(T &os, const RDPF &rdpf,
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
    os.write((const char *)&rdpf.unit_sum_inverse, sizeof(rdpf.unit_sum_inverse));
    os.write((const char *)&rdpf.scaled_sum, sizeof(rdpf.scaled_sum));
    os.write((const char *)&rdpf.scaled_xor, sizeof(rdpf.scaled_xor));

    return os;
}

// The ordinary << version never writes the expansion, since this is
// what we use to send DPFs over the network.
template <typename T>
T& operator<<(T &os, const RDPF &rdpf)
{
    return write_maybe_expanded(os, rdpf, false);
}

// I/O for RDPF Triples

// We never write RDPFTriples over the network, so always write
// the DPF expansions if they're available.
template <typename T>
T& operator<<(T &os, const RDPFTriple &rdpftrip)
{
    write_maybe_expanded(os, rdpftrip.dpf[0], true);
    write_maybe_expanded(os, rdpftrip.dpf[1], true);
    write_maybe_expanded(os, rdpftrip.dpf[2], true);
    nbits_t depth = rdpftrip.dpf[0].depth();
    os.write((const char *)&rdpftrip.as_target.ashare, BITBYTES(depth));
    os.write((const char *)&rdpftrip.xs_target.xshare, BITBYTES(depth));
    return os;
}

template <typename T>
T& operator>>(T &is, RDPFTriple &rdpftrip)
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
template <typename T>
T& operator<<(T &os, const RDPFPair &rdpfpair)
{
    write_maybe_expanded(os, rdpfpair.dpf[0], true);
    write_maybe_expanded(os, rdpfpair.dpf[1], true);
    return os;
}

template <typename T>
T& operator>>(T &is, RDPFPair &rdpfpair)
{
    is >> rdpfpair.dpf[0] >> rdpfpair.dpf[1];
    return is;
}
