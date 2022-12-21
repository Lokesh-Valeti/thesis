#include <bsd/stdlib.h> // arc4random_buf

#include "online.hpp"
#include "mpcops.hpp"
#include "rdpf.hpp"
#include "duoram.hpp"


static void online_test(MPCIO &mpcio, yield_t &yield,
    const PRACOptions &opts, char **args)
{
    nbits_t nbits = VALUE_BITS;

    if (*args) {
        nbits = atoi(*args);
    }

    size_t memsize = 9;

    MPCTIO tio(mpcio, 0);
    bool is_server = (mpcio.player == 2);

    RegAS *A = new RegAS[memsize];
    value_t V;
    RegBS F0, F1;
    RegXS X;

    if (!is_server) {
        A[0].randomize();
        A[1].randomize();
        F0.randomize();
        A[4].randomize();
        F1.randomize();
        A[6].randomize();
        A[7].randomize();
        X.randomize();
        arc4random_buf(&V, sizeof(V));
        printf("A:\n"); for (size_t i=0; i<memsize; ++i) printf("%3lu: %016lX\n", i, A[i].ashare);
        printf("V  : %016lX\n", V);
        printf("F0 : %01X\n", F0.bshare);
        printf("F1 : %01X\n", F1.bshare);
        printf("X  : %016lX\n", X.xshare);
    }
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_mul(tio, yield, A[2], A[0], A[1], nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_valuemul(tio, yield, A[3], V, nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_flagmult(tio, yield, A[5], F0, A[4], nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_oswap(tio, yield, A[6], A[7], F1, nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_xs_to_as(tio, yield, A[8], X, nbits);
        });
    run_coroutines(yield, coroutines);
    if (!is_server) {
        printf("\n");
        printf("A:\n"); for (size_t i=0; i<memsize; ++i) printf("%3lu: %016lX\n", i, A[i].ashare);
    }

    // Check the answers
    if (mpcio.player == 1) {
        tio.queue_peer(A, memsize*sizeof(RegAS));
        tio.queue_peer(&V, sizeof(V));
        tio.queue_peer(&F0, sizeof(RegBS));
        tio.queue_peer(&F1, sizeof(RegBS));
        tio.queue_peer(&X, sizeof(RegXS));
        tio.send();
    } else if (mpcio.player == 0) {
        RegAS *B = new RegAS[memsize];
        RegBS BF0, BF1;
        RegXS BX;
        value_t BV;
        value_t *S = new value_t[memsize];
        bit_t SF0, SF1;
        value_t SX;
        tio.recv_peer(B, memsize*sizeof(RegAS));
        tio.recv_peer(&BV, sizeof(BV));
        tio.recv_peer(&BF0, sizeof(RegBS));
        tio.recv_peer(&BF1, sizeof(RegBS));
        tio.recv_peer(&BX, sizeof(RegXS));
        for(size_t i=0; i<memsize; ++i) S[i] = A[i].ashare+B[i].ashare;
        SF0 = F0.bshare ^ BF0.bshare;
        SF1 = F1.bshare ^ BF1.bshare;
        SX = X.xshare ^ BX.xshare;
        printf("S:\n"); for (size_t i=0; i<memsize; ++i) printf("%3lu: %016lX\n", i, S[i]);
        printf("SF0: %01X\n", SF0);
        printf("SF1: %01X\n", SF1);
        printf("SX : %016lX\n", SX);
        printf("\n%016lx\n", S[0]*S[1]-S[2]);
        printf("%016lx\n", (V*BV)-S[3]);
        printf("%016lx\n", (SF0*S[4])-S[5]);
        printf("%016lx\n", S[8]-SX);
        delete[] B;
        delete[] S;
    }

    delete[] A;
}

static void lamport_test(MPCIO &mpcio, yield_t &yield,
    const PRACOptions &opts, char **args)
{
    // Create a bunch of threads and send a bunch of data to the other
    // peer, and receive their data.  If an arg is specified, repeat
    // that many times.  The Lamport clock at the end should be just the
    // number of repetitions.  Subsequent args are the chunk size and
    // the number of chunks per message
    size_t niters = 1;
    size_t chunksize = 1<<20;
    size_t numchunks = 1;
    if (*args) {
        niters = atoi(*args);
        ++args;
    }
    if (*args) {
        chunksize = atoi(*args);
        ++args;
    }
    if (*args) {
        numchunks = atoi(*args);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, niters, chunksize, numchunks] {
            MPCTIO tio(mpcio, thread_num);
            char *sendbuf = new char[chunksize];
            char *recvbuf = new char[chunksize*numchunks];
            for (size_t i=0; i<niters; ++i) {
                for (size_t chunk=0; chunk<numchunks; ++chunk) {
                    arc4random_buf(sendbuf, chunksize);
                    tio.queue_peer(sendbuf, chunksize);
                }
                tio.send();
                tio.recv_peer(recvbuf, chunksize*numchunks);
            }
            delete[] recvbuf;
            delete[] sendbuf;
        });
    }
    pool.join();
}

static void rdpf_test(MPCIO &mpcio, yield_t &yield,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, depth] {
            MPCTIO tio(mpcio, thread_num);
            size_t &op_counter = tio.aes_ops();
            if (mpcio.player == 2) {
                RDPFPair dp = tio.rdpfpair(depth);
                for (int i=0;i<2;++i) {
                    const RDPF &dpf = dp.dpf[i];
                    for (address_t x=0;x<(address_t(1)<<depth);++x) {
                        DPFnode leaf = dpf.leaf(x, op_counter);
                        RegBS ub = dpf.unit_bs(leaf);
                        RegAS ua = dpf.unit_as(leaf);
                        RegXS sx = dpf.scaled_xs(leaf);
                        RegAS sa = dpf.scaled_as(leaf);
                        printf("%04x %x %016lx %016lx %016lx\n", x,
                            ub.bshare, ua.ashare, sx.xshare, sa.ashare);
                    }
                    printf("\n");
                }
            } else {
                RDPFTriple dt = tio.rdpftriple(depth);
                for (int i=0;i<3;++i) {
                    const RDPF &dpf = dt.dpf[i];
                    RegXS peer_scaled_xor;
                    RegAS peer_scaled_sum;
                    if (mpcio.player == 1) {
                        tio.iostream_peer() << dpf.scaled_xor << dpf.scaled_sum;
                    } else {
                        tio.iostream_peer() >> peer_scaled_xor >> peer_scaled_sum;
                        peer_scaled_sum += dpf.scaled_sum;
                        peer_scaled_xor ^= dpf.scaled_xor;
                    }
                    for (address_t x=0;x<(address_t(1)<<depth);++x) {
                        DPFnode leaf = dpf.leaf(x, op_counter);
                        RegBS ub = dpf.unit_bs(leaf);
                        RegAS ua = dpf.unit_as(leaf);
                        RegXS sx = dpf.scaled_xs(leaf);
                        RegAS sa = dpf.scaled_as(leaf);
                        printf("%04x %x %016lx %016lx %016lx\n", x,
                            ub.bshare, ua.ashare, sx.xshare, sa.ashare);
                        if (mpcio.player == 1) {
                            tio.iostream_peer() << ub << ua << sx << sa;
                        } else {
                            RegBS peer_ub;
                            RegAS peer_ua;
                            RegXS peer_sx;
                            RegAS peer_sa;
                            tio.iostream_peer() >> peer_ub >> peer_ua >>
                                peer_sx >> peer_sa;
                            ub ^= peer_ub;
                            ua += peer_ua;
                            sx ^= peer_sx;
                            sa += peer_sa;
                            if (ub.bshare || ua.ashare || sx.xshare || sa.ashare) {
                                printf("**** %x %016lx %016lx %016lx\n",
                                    ub.bshare, ua.ashare, sx.xshare, sa.ashare);
                                printf("SCALE                   %016lx %016lx\n",
                                    peer_scaled_xor.xshare, peer_scaled_sum.ashare);
                            }
                        }
                    }
                    printf("\n");
                }
            }
            tio.send();
        });
    }
    pool.join();
}

static void rdpf_timing(MPCIO &mpcio, yield_t &yield,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, depth] {
            MPCTIO tio(mpcio, thread_num);
            size_t &op_counter = tio.aes_ops();
            if (mpcio.player == 2) {
                RDPFPair dp = tio.rdpfpair(depth);
                for (int i=0;i<2;++i) {
                    RDPF &dpf = dp.dpf[i];
                    dpf.expand(op_counter);
                    RegXS scaled_xor;
                    for (address_t x=0;x<(address_t(1)<<depth);++x) {
                        DPFnode leaf = dpf.leaf(x, op_counter);
                        RegXS sx = dpf.scaled_xs(leaf);
                        scaled_xor ^= sx;
                    }
                    printf("%016lx\n%016lx\n", scaled_xor.xshare,
                        dpf.scaled_xor.xshare);
                    printf("\n");
                }
            } else {
                RDPFTriple dt = tio.rdpftriple(depth);
                for (int i=0;i<3;++i) {
                    RDPF &dpf = dt.dpf[i];
                    dpf.expand(op_counter);
                    RegXS scaled_xor;
                    for (address_t x=0;x<(address_t(1)<<depth);++x) {
                        DPFnode leaf = dpf.leaf(x, op_counter);
                        RegXS sx = dpf.scaled_xs(leaf);
                        scaled_xor ^= sx;
                    }
                    printf("%016lx\n%016lx\n", scaled_xor.xshare,
                        dpf.scaled_xor.xshare);
                    printf("\n");
                }
            }
            tio.send();
        });
    }
    pool.join();
}

static void rdpfeval_timing(MPCIO &mpcio, yield_t &yield,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;
    address_t start=0;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    if (*args) {
        start = atoi(*args);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, depth, start] {
            MPCTIO tio(mpcio, thread_num);
            size_t &op_counter = tio.aes_ops();
            if (mpcio.player == 2) {
                RDPFPair dp = tio.rdpfpair(depth);
                for (int i=0;i<2;++i) {
                    RDPF &dpf = dp.dpf[i];
                    RegXS scaled_xor;
                    auto ev = StreamEval(dpf, start, op_counter, false);
                    for (address_t x=0;x<(address_t(1)<<depth);++x) {
                        DPFnode leaf = ev.next();
                        RegXS sx = dpf.scaled_xs(leaf);
                        scaled_xor ^= sx;
                    }
                    printf("%016lx\n%016lx\n", scaled_xor.xshare,
                        dpf.scaled_xor.xshare);
                    printf("\n");
                }
            } else {
                RDPFTriple dt = tio.rdpftriple(depth);
                for (int i=0;i<3;++i) {
                    RDPF &dpf = dt.dpf[i];
                    RegXS scaled_xor;
                    auto ev = StreamEval(dpf, start, op_counter, false);
                    for (address_t x=0;x<(address_t(1)<<depth);++x) {
                        DPFnode leaf = ev.next();
                        RegXS sx = dpf.scaled_xs(leaf);
                        scaled_xor ^= sx;
                    }
                    printf("%016lx\n%016lx\n", scaled_xor.xshare,
                        dpf.scaled_xor.xshare);
                    printf("\n");
                }
            }
            tio.send();
        });
    }
    pool.join();
}

static void tupleeval_timing(MPCIO &mpcio, yield_t &yield,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;
    address_t start=0;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    if (*args) {
        start = atoi(*args);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, depth, start] {
            MPCTIO tio(mpcio, thread_num);
            size_t &op_counter = tio.aes_ops();
            if (mpcio.player == 2) {
                RDPFPair dp = tio.rdpfpair(depth);
                RegXS scaled_xor0, scaled_xor1;
                auto ev = StreamEval(dp, start, op_counter, false);
                for (address_t x=0;x<(address_t(1)<<depth);++x) {
                    auto [L0, L1] = ev.next();
                    RegXS sx0 = dp.dpf[0].scaled_xs(L0);
                    RegXS sx1 = dp.dpf[1].scaled_xs(L1);
                    scaled_xor0 ^= sx0;
                    scaled_xor1 ^= sx1;
                }
                printf("%016lx\n%016lx\n", scaled_xor0.xshare,
                    dp.dpf[0].scaled_xor.xshare);
                printf("\n");
                printf("%016lx\n%016lx\n", scaled_xor1.xshare,
                    dp.dpf[1].scaled_xor.xshare);
                printf("\n");
            } else {
                RDPFTriple dt = tio.rdpftriple(depth);
                RegXS scaled_xor0, scaled_xor1, scaled_xor2;
                auto ev = StreamEval(dt, start, op_counter, false);
                for (address_t x=0;x<(address_t(1)<<depth);++x) {
                    auto [L0, L1, L2] = ev.next();
                    RegXS sx0 = dt.dpf[0].scaled_xs(L0);
                    RegXS sx1 = dt.dpf[1].scaled_xs(L1);
                    RegXS sx2 = dt.dpf[2].scaled_xs(L2);
                    scaled_xor0 ^= sx0;
                    scaled_xor1 ^= sx1;
                    scaled_xor2 ^= sx2;
                }
                printf("%016lx\n%016lx\n", scaled_xor0.xshare,
                    dt.dpf[0].scaled_xor.xshare);
                printf("\n");
                printf("%016lx\n%016lx\n", scaled_xor1.xshare,
                    dt.dpf[1].scaled_xor.xshare);
                printf("\n");
                printf("%016lx\n%016lx\n", scaled_xor2.xshare,
                    dt.dpf[2].scaled_xor.xshare);
                printf("\n");
            }
            tio.send();
        });
    }
    pool.join();
}

static void duoram_test(MPCIO &mpcio, yield_t &yield,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, &yield, thread_num, depth] {
            MPCTIO tio(mpcio, thread_num);
            // size_t &op_counter = tio.aes_ops();
            Duoram<RegAS> oram(mpcio.player, size_t(1)<<depth);
            printf("%ld\n", oram.size());
            auto A = oram.flat(tio, yield);
            RegAS aidx;
            aidx.randomize(depth);
            RegXS xidx;
            xidx.randomize(depth);
            size_t eidx = arc4random();
            eidx &= (size_t(1)<<depth)-1;
            RegAS Aa = A[aidx];
            auto Ax = A[xidx];
            auto Ae = A[eidx];
            printf("%ld %ld\n", A.size(), Aa.ashare);
            tio.send();
        });
    }
    pool.join();
}

void online_main(MPCIO &mpcio, const PRACOptions &opts, char **args)
{
    // Run everything inside a coroutine so that simple tests don't have
    // to start one themselves
    MPCTIO tio(mpcio, 0);
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(
        [&](yield_t &yield) {
            if (!*args) {
                std::cerr << "Mode is required as the first argument when not preprocessing.\n";
                return;
            } else if (!strcmp(*args, "test")) {
                ++args;
                online_test(mpcio, yield, opts, args);
            } else if (!strcmp(*args, "lamporttest")) {
                ++args;
                lamport_test(mpcio, yield, opts, args);
            } else if (!strcmp(*args, "rdpftest")) {
                ++args;
                rdpf_test(mpcio, yield, opts, args);
            } else if (!strcmp(*args, "rdpftime")) {
                ++args;
                rdpf_timing(mpcio, yield, opts, args);
            } else if (!strcmp(*args, "evaltime")) {
                ++args;
                rdpfeval_timing(mpcio, yield, opts, args);
            } else if (!strcmp(*args, "tupletime")) {
                ++args;
                tupleeval_timing(mpcio, yield, opts, args);
            } else if (!strcmp(*args, "duotest")) {
                ++args;
                duoram_test(mpcio, yield, opts, args);
            } else {
                std::cerr << "Unknown mode " << *args << "\n";
            }
        });
    run_coroutines(tio, coroutines);
}
