#ifndef __OPTIONS_HPP__
#define __OPTIONS_HPP__

#include "types.hpp"

struct PRACOptions {
    ProcessingMode mode;
    int num_threads;
    bool expand_rdpfs;
    bool use_xor_db;
    bool append_to_files;

    PRACOptions() : mode(MODE_ONLINE), num_threads(1),
        expand_rdpfs(false), use_xor_db(false), append_to_files(false) {}
};

#endif
