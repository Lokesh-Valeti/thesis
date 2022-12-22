#ifndef __OPTIONS_HPP__
#define __OPTIONS_HPP__

struct PRACOptions {
    bool preprocessing;
    int num_threads;
    bool expand_rdpfs;
    bool use_xor_db;

    PRACOptions() : preprocessing(false), num_threads(1),
        expand_rdpfs(true), use_xor_db(false) {}
};

#endif
