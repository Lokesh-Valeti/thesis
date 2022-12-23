// Templated method implementations for cdpf.hpp

// I/O for CDPFs

// Read the DPF from the output stream.  You can use this to read DPFs
// from files or from the network.
template <typename T>
T& operator>>(T &is, CDPF &cdpf)
{
    is.read((char *)&cdpf.seed, sizeof(cdpf.seed));
    cdpf.whichhalf = get_lsb(cdpf.seed);
    uint8_t depth = VALUE_BITS - 7;
    cdpf.cw.clear();
    for (uint8_t i=0; i<depth; ++i) {
        DPFnode cw;
        is.read((char *)&cw, sizeof(cw));
        cdpf.cw.push_back(cw);
    }
    value_t cfbits = 0;
    is.read((char *)&cfbits, BITBYTES(depth));
    cdpf.cfbits = cfbits;
    is.read((char *)&cdpf.leaf_cwr, sizeof(cdpf.leaf_cwr));
    is.read((char *)&cdpf.as_target, sizeof(cdpf.as_target));
    is.read((char *)&cdpf.xs_target, sizeof(cdpf.xs_target));

    return is;
}

// Write the DPF to the output stream.  You can use this to write DPFs
// to files or over the network.
template <typename T>
T& operator<<(T &os, const CDPF &cdpf)
{
    os.write((const char *)&cdpf.seed, sizeof(cdpf.seed));
    uint8_t depth = VALUE_BITS - 7;
    for (uint8_t i=0; i<depth; ++i) {
        os.write((const char *)&cdpf.cw[i], sizeof(cdpf.cw[i]));
    }
    os.write((const char *)&cdpf.cfbits, BITBYTES(depth));
    os.write((const char *)&cdpf.leaf_cwr, sizeof(cdpf.leaf_cwr));
    os.write((const char *)&cdpf.as_target, sizeof(cdpf.as_target));
    os.write((const char *)&cdpf.xs_target, sizeof(cdpf.xs_target));

    return os;
}
