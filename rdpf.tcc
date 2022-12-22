// Templated method implementations for rdpf.hpp

// Create a StreamEval object that will start its output at index start.
// It will wrap around to 0 when it hits 2^depth.  If use_expansion
// is true, then if the DPF has been expanded, just output values
// from that.  If use_expansion=false or if the DPF has not been
// expanded, compute the values on the fly.
template <typename T>
StreamEval<T>::StreamEval(const T &rdpf, address_t start,
    size_t &op_counter, bool use_expansion) : rdpf(rdpf),
    op_counter(op_counter), use_expansion(use_expansion)
{
    depth = rdpf.depth();
    // Prevent overflow of 1<<depth
    if (depth < ADDRESS_MAX_BITS) {
        indexmask = (address_t(1)<<depth)-1;
    } else {
        indexmask = ~0;
    }
    start &= indexmask;
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
        path[i] = rdpf.descend(path[i-1], i-1, dir, op_counter);
    }
}

template <typename T>
typename T::node StreamEval<T>::next()
{
    if (use_expansion && rdpf.has_expansion()) {
        // Just use the precomputed values
        typename T::node leaf = rdpf.get_expansion(nextindex);
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
            nextindex & (address_t(1) << how_many_1_bits);
        path[depth-how_many_1_bits] =
            rdpf.descend(path[depth-how_many_1_bits-1],
                depth-how_many_1_bits-1, top_changed_bit, op_counter);
        for (nbits_t i = depth-how_many_1_bits; i < depth-1; ++i) {
            path[i+1] = rdpf.descend(path[i], i, 0, op_counter);
        }
    }
    typename T::node leaf = rdpf.descend(path[depth-1], depth-1,
        nextindex & 1, op_counter);
    pathindex = nextindex;
    nextindex = (nextindex + 1) & indexmask;
    return leaf;
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
    uint8_t depth;
    // The whichhalf bit is the high bit of depth
    is.read((char *)&depth, sizeof(depth));
    rdpf.whichhalf = !!(depth & 0x80);
    depth &= 0x7f;
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
    // The whichhalf bit is the high bit of depth
    // If we're writing an expansion, add 64 to depth as well
    uint8_t whichhalf_and_depth = depth |
        (uint8_t(rdpf.whichhalf)<<7);
    bool write_expansion = false;
    if (expanded && rdpf.expansion.size() == (size_t(1)<<depth)) {
        write_expansion = true;
        whichhalf_and_depth += 64;
    }
    os.write((const char *)&whichhalf_and_depth,
        sizeof(whichhalf_and_depth));
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
