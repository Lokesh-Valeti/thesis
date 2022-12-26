#ifndef __OPTIONS_HPP__
#define __OPTIONS_HPP__

#include "types.hpp"

struct PRACOptions {
    ProcessingMode mode;
    int num_threads;
    bool expand_rdpfs;
    bool use_xor_db;

    PRACOptions() : mode(MODE_ONLINE), num_threads(1),
        expand_rdpfs(true), use_xor_db(false) {}
};

#endif
