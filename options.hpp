#ifndef __OPTIONS_HPP__
#define __OPTIONS_HPP__

struct PRACOptions {
    bool preprocessing;
    int num_threads;
    bool expand_rdpfs;

    PRACOptions() : preprocessing(false), num_threads(1),
        expand_rdpfs(true) {}
};

#endif
