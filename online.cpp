#include <bsd/stdlib.h> // arc4random_buf

#include "online.hpp"
#include "mpcops.hpp"
#include "rdpf.hpp"
#include "duoram.hpp"
#include "cdpf.hpp"


static void online_test(MPCIO &mpcio,
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
        [&tio, &A, nbits](yield_t &yield) {
            mpc_mul(tio, yield, A[2], A[0], A[1], nbits);
        });
    coroutines.emplace_back(
        [&tio, &A, V, nbits](yield_t &yield) {
            mpc_valuemul(tio, yield, A[3], V, nbits);
        });
    coroutines.emplace_back(
        [&tio, &A, &F0, nbits](yield_t &yield) {
            mpc_flagmult(tio, yield, A[5], F0, A[4], nbits);
        });
    coroutines.emplace_back(
        [&tio, &A, &F1, nbits](yield_t &yield) {
            mpc_oswap(tio, yield, A[6], A[7], F1, nbits);
        });
    coroutines.emplace_back(
        [&tio, &A, &X, nbits](yield_t &yield) {
            mpc_xs_to_as(tio, yield, A[8], X, nbits);
        });
    run_coroutines(tio, coroutines);
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

static void lamport_test(MPCIO &mpcio,
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

static void rdpf_test(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;
    size_t num_iters = 1;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    if (*args) {
        num_iters = atoi(*args);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, depth, num_iters] {
            MPCTIO tio(mpcio, thread_num);
            run_coroutines(tio, [&tio, depth, num_iters] (yield_t &yield) {
                size_t &aes_ops = tio.aes_ops();
                for (size_t iter=0; iter < num_iters; ++iter) {
                    if (tio.player() == 2) {
                        RDPFPair dp = tio.rdpfpair(yield, depth);
                        for (int i=0;i<2;++i) {
                            const RDPF &dpf = dp.dpf[i];
                            for (address_t x=0;x<(address_t(1)<<depth);++x) {
                                DPFnode leaf = dpf.leaf(x, aes_ops);
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
                        RDPFTriple dt = tio.rdpftriple(yield, depth);
                        for (int i=0;i<3;++i) {
                            const RDPF &dpf = dt.dpf[i];
                            RegXS peer_scaled_xor;
                            RegAS peer_scaled_sum;
                            if (tio.player() == 1) {
                                tio.iostream_peer() << dpf.scaled_xor << dpf.scaled_sum;
                            } else {
                                tio.iostream_peer() >> peer_scaled_xor >> peer_scaled_sum;
                                peer_scaled_sum += dpf.scaled_sum;
                                peer_scaled_xor ^= dpf.scaled_xor;
                            }
                            for (address_t x=0;x<(address_t(1)<<depth);++x) {
                                DPFnode leaf = dpf.leaf(x, aes_ops);
                                RegBS ub = dpf.unit_bs(leaf);
                                RegAS ua = dpf.unit_as(leaf);
                                RegXS sx = dpf.scaled_xs(leaf);
                                RegAS sa = dpf.scaled_as(leaf);
                                printf("%04x %x %016lx %016lx %016lx\n", x,
                                    ub.bshare, ua.ashare, sx.xshare, sa.ashare);
                                if (tio.player() == 1) {
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
                }
            });
        });
    }
    pool.join();
}

static void rdpf_timing(MPCIO &mpcio,
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
            run_coroutines(tio, [&tio, depth] (yield_t &yield) {
                size_t &aes_ops = tio.aes_ops();
                if (tio.player() == 2) {
                    RDPFPair dp = tio.rdpfpair(yield, depth);
                    for (int i=0;i<2;++i) {
                        RDPF &dpf = dp.dpf[i];
                        dpf.expand(aes_ops);
                        RegXS scaled_xor;
                        for (address_t x=0;x<(address_t(1)<<depth);++x) {
                            DPFnode leaf = dpf.leaf(x, aes_ops);
                            RegXS sx = dpf.scaled_xs(leaf);
                            scaled_xor ^= sx;
                        }
                        printf("%016lx\n%016lx\n", scaled_xor.xshare,
                            dpf.scaled_xor.xshare);
                        printf("\n");
                    }
                } else {
                    RDPFTriple dt = tio.rdpftriple(yield, depth);
                    for (int i=0;i<3;++i) {
                        RDPF &dpf = dt.dpf[i];
                        dpf.expand(aes_ops);
                        RegXS scaled_xor;
                        for (address_t x=0;x<(address_t(1)<<depth);++x) {
                            DPFnode leaf = dpf.leaf(x, aes_ops);
                            RegXS sx = dpf.scaled_xs(leaf);
                            scaled_xor ^= sx;
                        }
                        printf("%016lx\n%016lx\n", scaled_xor.xshare,
                            dpf.scaled_xor.xshare);
                        printf("\n");
                    }
                }
            });
        });
    }
    pool.join();
}

static value_t parallel_streameval_rdpf(MPCIO &mpcio, const RDPF &dpf,
    address_t start, int num_threads)
{
    RegXS scaled_xor[num_threads];
    boost::asio::thread_pool pool(num_threads);
    address_t totsize = (address_t(1)<<dpf.depth());
    address_t threadstart = start;
    address_t threadchunk = totsize / num_threads;
    address_t threadextra = totsize % num_threads;
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        address_t threadsize = threadchunk + (address_t(thread_num) < threadextra);
        boost::asio::post(pool,
            [&mpcio, &dpf, &scaled_xor, thread_num, threadstart, threadsize] {
                MPCTIO tio(mpcio, thread_num);
//printf("Thread %d from %X for %X\n", thread_num, threadstart, threadsize);
                RegXS local_xor;
                size_t local_aes_ops = 0;
                auto ev = StreamEval(dpf, threadstart, 0, local_aes_ops);
                for (address_t x=0;x<threadsize;++x) {
//if (x%0x10000 == 0) printf("%d", thread_num);
                    DPFnode leaf = ev.next();
                    local_xor ^= dpf.scaled_xs(leaf);
                }
                scaled_xor[thread_num] = local_xor;
                tio.aes_ops() += local_aes_ops;
//printf("Thread %d complete\n", thread_num);
            });
        threadstart = (threadstart + threadsize) % totsize;
    }
    pool.join();
    RegXS res;
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        res ^= scaled_xor[thread_num];
    }
    return res.xshare;
}

static void rdpfeval_timing(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;
    address_t start=0;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    if (*args) {
        start = strtoull(*args, NULL, 16);
        ++args;
    }

    int num_threads = opts.num_threads;
    MPCTIO tio(mpcio, 0, num_threads);
    run_coroutines(tio, [&mpcio, &tio, depth, start, num_threads] (yield_t &yield) {
        if (tio.player() == 2) {
            RDPFPair dp = tio.rdpfpair(yield, depth);
            for (int i=0;i<2;++i) {
                RDPF &dpf = dp.dpf[i];
                value_t scaled_xor =
                    parallel_streameval_rdpf(mpcio, dpf, start, num_threads);
                printf("%016lx\n%016lx\n", scaled_xor,
                    dpf.scaled_xor.xshare);
                printf("\n");
            }
        } else {
            RDPFTriple dt = tio.rdpftriple(yield, depth);
            for (int i=0;i<3;++i) {
                RDPF &dpf = dt.dpf[i];
                value_t scaled_xor =
                    parallel_streameval_rdpf(mpcio, dpf, start, num_threads);
                printf("%016lx\n%016lx\n", scaled_xor,
                    dpf.scaled_xor.xshare);
                printf("\n");
            }
        }
    });
}

static void par_rdpfeval_timing(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;
    address_t start=0;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    if (*args) {
        start = strtoull(*args, NULL, 16);
        ++args;
    }

    int num_threads = opts.num_threads;
    MPCTIO tio(mpcio, 0, num_threads);
    run_coroutines(tio, [&tio, depth, start, num_threads] (yield_t &yield) {
        if (tio.player() == 2) {
            RDPFPair dp = tio.rdpfpair(yield, depth);
            for (int i=0;i<2;++i) {
                RDPF &dpf = dp.dpf[i];
                nbits_t depth = dpf.depth();
                auto pe = ParallelEval(dpf, start, 0,
                    address_t(1)<<depth, num_threads, tio.aes_ops());
                RegXS result, init;
                result = pe.reduce(init, [&dpf] (int thread_num,
                        address_t i, const RDPF::node &leaf) {
                    return dpf.scaled_xs(leaf);
                });
                printf("%016lx\n%016lx\n", result.xshare,
                    dpf.scaled_xor.xshare);
                printf("\n");
            }
        } else {
            RDPFTriple dt = tio.rdpftriple(yield, depth);
            for (int i=0;i<3;++i) {
                RDPF &dpf = dt.dpf[i];
                nbits_t depth = dpf.depth();
                auto pe = ParallelEval(dpf, start, 0,
                    address_t(1)<<depth, num_threads, tio.aes_ops());
                RegXS result, init;
                result = pe.reduce(init, [&dpf] (int thread_num,
                        address_t i, const RDPF::node &leaf) {
                    return dpf.scaled_xs(leaf);
                });
                printf("%016lx\n%016lx\n", result.xshare,
                    dpf.scaled_xor.xshare);
                printf("\n");
            }
        }
    });
}

static void tupleeval_timing(MPCIO &mpcio,
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
    MPCTIO tio(mpcio, 0, num_threads);
    run_coroutines(tio, [&tio, depth, start] (yield_t &yield) {
        size_t &aes_ops = tio.aes_ops();
        if (tio.player() == 2) {
            RDPFPair dp = tio.rdpfpair(yield, depth);
            RegXS scaled_xor0, scaled_xor1;
            auto ev = StreamEval(dp, start, 0, aes_ops, false);
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
            RDPFTriple dt = tio.rdpftriple(yield, depth);
            RegXS scaled_xor0, scaled_xor1, scaled_xor2;
            auto ev = StreamEval(dt, start, 0, aes_ops, false);
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
    });
}

static void par_tupleeval_timing(MPCIO &mpcio,
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
    MPCTIO tio(mpcio, 0, num_threads);
    run_coroutines(tio, [&tio, depth, start, num_threads] (yield_t &yield) {
        size_t &aes_ops = tio.aes_ops();
        if (tio.player() == 2) {
            RDPFPair dp = tio.rdpfpair(yield, depth);
            auto pe = ParallelEval(dp, start, 0, address_t(1)<<depth,
                num_threads, aes_ops);
            using V = std::tuple<RegXS,RegXS>;
            V result, init;
            result = pe.reduce(init, [&dp] (int thread_num, address_t i,
                    const RDPFPair::node &leaf) {
                return dp.scaled<RegXS>(leaf);
            });
            printf("%016lx\n%016lx\n", std::get<0>(result).xshare,
                dp.dpf[0].scaled_xor.xshare);
            printf("\n");
            printf("%016lx\n%016lx\n", std::get<1>(result).xshare,
                dp.dpf[1].scaled_xor.xshare);
            printf("\n");
        } else {
            RDPFTriple dt = tio.rdpftriple(yield, depth);
            auto pe = ParallelEval(dt, start, 0, address_t(1)<<depth,
                num_threads, aes_ops);
            using V = std::tuple<RegXS,RegXS,RegXS>;
            V result, init;
            result = pe.reduce(init, [&dt] (int thread_num, address_t i,
                    const RDPFTriple::node &leaf) {
                return dt.scaled<RegXS>(leaf);
            });
            printf("%016lx\n%016lx\n", std::get<0>(result).xshare,
                dt.dpf[0].scaled_xor.xshare);
            printf("\n");
            printf("%016lx\n%016lx\n", std::get<1>(result).xshare,
                dt.dpf[1].scaled_xor.xshare);
            printf("\n");
            printf("%016lx\n%016lx\n", std::get<2>(result).xshare,
                dt.dpf[2].scaled_xor.xshare);
            printf("\n");
        }
    });
}

// T is RegAS or RegXS for additive or XOR shared database respectively
template <typename T>
static void duoram_test(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=6;
    address_t share=arc4random();

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    if (*args) {
        share = atoi(*args);
        ++args;
    }
    share &= ((address_t(1)<<depth)-1);

    MPCTIO tio(mpcio, 0, opts.num_threads);
    run_coroutines(tio, [&tio, depth, share] (yield_t &yield) {
        size_t size = size_t(1)<<depth;
        // size_t &aes_ops = tio.aes_ops();
        Duoram<T> oram(tio.player(), size);
        auto A = oram.flat(tio, yield);
        RegAS aidx, aidx2, aidx3;
        aidx.ashare = share;
        aidx2.ashare = share + tio.player();
        aidx3.ashare = share + 1;
        T M;
        if (tio.player() == 0) {
            M.set(0xbabb0000);
        } else {
            M.set(0x0000a66e);
        }
        RegXS xidx;
        xidx.xshare = share;
        T N;
        if (tio.player() == 0) {
            N.set(0xdead0000);
        } else {
            N.set(0x0000beef);
        }
        // Writing and reading with additively shared indices
        printf("Updating\n");
        A[aidx] += M;
        printf("Reading\n");
        T Aa = A[aidx];
        // Writing and reading with XOR shared indices
        printf("Updating\n");
        A[xidx] += N;
        printf("Reading\n");
        T Ax = A[xidx];
        T Ae;
        // Writing and reading with explicit indices
        if (depth > 2) {
            A[5] += Aa;
            Ae = A[6];
        }

        // Simultaneous independent reads
        std::vector<T> Av;
        Av.resize(3);
        std::vector<coro_t> coroutines;
        run_coroutines(yield,
            [&A, &Av, &aidx] (yield_t &yield) {
                auto Acoro = A.context(yield);
                Av[0] = Acoro[aidx];
            },
            [&A, &Av, &aidx2] (yield_t &yield) {
                auto Acoro = A.context(yield);
                Av[1] = Acoro[aidx2];
            },
            [&A, &Av, &aidx3] (yield_t &yield) {
                auto Acoro = A.context(yield);
                Av[2] = Acoro[aidx3];
            });

        if (depth <= 10) {
            oram.dump();
            auto check = A.reconstruct();
            if (tio.player() == 0) {
                for (address_t i=0;i<size;++i) {
                    printf("%04x %016lx\n", i, check[i].share());
                }
            }
        }
        auto checkread = A.reconstruct(Aa);
        auto checkreade = A.reconstruct(Ae);
        auto checkreadx = A.reconstruct(Ax);
        if (tio.player() == 0) {
            printf("Read AS value = %016lx\n", checkread.share());
            printf("Read AX value = %016lx\n", checkreadx.share());
            printf("Read Ex value = %016lx\n", checkreade.share());
        }
        for (auto &v : Av) {
            auto checkv = A.reconstruct(v);
            if (tio.player() == 0) {
                printf("Read Av value = %016lx\n", checkv.share());
            }
        }
    });
}

static void cdpf_test(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    value_t query, target;
    int iters = 1;
    arc4random_buf(&query, sizeof(query));
    arc4random_buf(&target, sizeof(target));

    if (*args) {
        query = strtoull(*args, NULL, 16);
        ++args;
    }
    if (*args) {
        target = strtoull(*args, NULL, 16);
        ++args;
    }
    if (*args) {
        iters = atoi(*args);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, query, target, iters] {
            MPCTIO tio(mpcio, thread_num);
            run_coroutines(tio, [&tio, query, target, iters] (yield_t &yield) {
                size_t &aes_ops = tio.aes_ops();
                for (int i=0;i<iters;++i) {
                    if (tio.player() == 2) {
                        tio.cdpf(yield);
                        auto [ dpf0, dpf1 ] = CDPF::generate(target, aes_ops);
                        DPFnode leaf0 = dpf0.leaf(query, aes_ops);
                        DPFnode leaf1 = dpf1.leaf(query, aes_ops);
                        printf("DPFXOR_{%016lx}(%016lx} = ", target, query);
                        dump_node(leaf0 ^ leaf1);
                    } else {
                        CDPF dpf = tio.cdpf(yield);
                        printf("ashare = %016lX\nxshare = %016lX\n",
                            dpf.as_target.ashare, dpf.xs_target.xshare);
                        DPFnode leaf = dpf.leaf(query, aes_ops);
                        printf("DPF(%016lx) = ", query);
                        dump_node(leaf);
                        if (tio.player() == 1) {
                            tio.iostream_peer() << leaf;
                        } else {
                            DPFnode peerleaf;
                            tio.iostream_peer() >> peerleaf;
                            printf("XOR = ");
                            dump_node(leaf ^ peerleaf);
                        }
                    }
                }
            });
        });
    }
    pool.join();
}

static int compare_test_one(MPCTIO &tio, yield_t &yield,
    value_t target, value_t x)
{
    int player = tio.player();
    size_t &aes_ops = tio.aes_ops();
    int res = 1;
    if (player == 2) {
        // Create a CDPF pair with the given target
        auto [dpf0, dpf1] = CDPF::generate(target, aes_ops);
        // Send it and a share of x to the computational parties
        RegAS x0, x1;
        x0.randomize();
        x1.set(x-x0.share());
        tio.iostream_p0() << dpf0 << x0;
        tio.iostream_p1() << dpf1 << x1;
    } else {
        CDPF dpf;
        RegAS xsh;
        tio.iostream_server() >> dpf >> xsh;
        auto [lt, eq, gt] = dpf.compare(tio, yield, xsh, aes_ops);
        printf("%016lx %016lx %d %d %d ", target, x, lt.bshare,
            eq.bshare, gt.bshare);
        // Check the answer
        if (player == 1) {
            tio.iostream_peer() << xsh << lt << eq << gt;
        } else {
            RegAS peer_xsh;
            RegBS peer_lt, peer_eq, peer_gt;
            tio.iostream_peer() >> peer_xsh >> peer_lt >> peer_eq >> peer_gt;
            lt ^= peer_lt;
            eq ^= peer_eq;
            gt ^= peer_gt;
            xsh += peer_xsh;
            int lti = int(lt.bshare);
            int eqi = int(eq.bshare);
            int gti = int(gt.bshare);
            x = xsh.share();
            printf(": %d %d %d ", lti, eqi, gti);
            bool signbit = (x >> 63);
            if (lti + eqi + gti != 1) {
                printf("INCONSISTENT");
                res = 0;
            } else if (x == 0 && eqi) {
                printf("=");
            } else if (!signbit && gti) {
                printf(">");
            } else if (signbit && lti) {
                printf("<");
            } else {
                printf("INCORRECT");
                res = 0;
            }
        }
        printf("\n");
    }
    return res;
}

static int compare_test_target(MPCTIO &tio, yield_t &yield,
    value_t target, value_t x)
{
    int res = 1;
    res &= compare_test_one(tio, yield, target, x);
    res &= compare_test_one(tio, yield, target, 0);
    res &= compare_test_one(tio, yield, target, 1);
    res &= compare_test_one(tio, yield, target, 15);
    res &= compare_test_one(tio, yield, target, 16);
    res &= compare_test_one(tio, yield, target, 17);
    res &= compare_test_one(tio, yield, target, -1);
    res &= compare_test_one(tio, yield, target, -15);
    res &= compare_test_one(tio, yield, target, -16);
    res &= compare_test_one(tio, yield, target, -17);
    res &= compare_test_one(tio, yield, target, (value_t(1)<<63));
    res &= compare_test_one(tio, yield, target, (value_t(1)<<63)+1);
    res &= compare_test_one(tio, yield, target, (value_t(1)<<63)-1);
    return res;
}

static void compare_test(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    value_t target, x;
    arc4random_buf(&target, sizeof(target));
    arc4random_buf(&x, sizeof(x));

    if (*args) {
        target = strtoull(*args, NULL, 16);
        ++args;
    }
    if (*args) {
        x = strtoull(*args, NULL, 16);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, target, x] {
            MPCTIO tio(mpcio, thread_num);
            run_coroutines(tio, [&tio, target, x] (yield_t &yield) {
                int res = 1;
                res &= compare_test_target(tio, yield, target, x);
                res &= compare_test_target(tio, yield, 0, x);
                res &= compare_test_target(tio, yield, 1, x);
                res &= compare_test_target(tio, yield, 15, x);
                res &= compare_test_target(tio, yield, 16, x);
                res &= compare_test_target(tio, yield, 17, x);
                res &= compare_test_target(tio, yield, -1, x);
                res &= compare_test_target(tio, yield, -15, x);
                res &= compare_test_target(tio, yield, -16, x);
                res &= compare_test_target(tio, yield, -17, x);
                res &= compare_test_target(tio, yield, (value_t(1)<<63), x);
                res &= compare_test_target(tio, yield, (value_t(1)<<63)+1, x);
                res &= compare_test_target(tio, yield, (value_t(1)<<63)-1, x);
                if (tio.player() == 0) {
                    if (res == 1) {
                        printf("All tests passed!\n");
                    } else {
                        printf("TEST FAILURES\n");
                    }
                }
            });
        });
    }
    pool.join();
}

static void sort_test(MPCIO &mpcio,
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
            run_coroutines(tio, [&tio, depth] (yield_t &yield) {
                address_t size = address_t(1)<<depth;
                // size_t &aes_ops = tio.aes_ops();
                Duoram<RegAS> oram(tio.player(), size);
                auto A = oram.flat(tio, yield);
                A.explicitonly(true);
                // Initialize the memory to random values in parallel
                std::vector<coro_t> coroutines;
                for (address_t i=0; i<size; ++i) {
                    coroutines.emplace_back(
                        [&A, i](yield_t &yield) {
                            auto Acoro = A.context(yield);
                            RegAS v;
                            v.randomize(62);
                            Acoro[i] += v;
                        });
                }
                run_coroutines(yield, coroutines);
                A.bitonic_sort(0, depth);
                if (depth <= 10) {
                    oram.dump();
                    auto check = A.reconstruct();
                    if (tio.player() == 0) {
                        for (address_t i=0;i<size;++i) {
                            printf("%04x %016lx\n", i, check[i].share());
                        }
                    }
                }
            });
        });
    }
    pool.join();
}

static void bsearch_test(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    value_t target;
    arc4random_buf(&target, sizeof(target));
    target >>= 1;
    nbits_t depth=6;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    if (*args) {
        target = strtoull(*args, NULL, 16);
        ++args;
    }

    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, depth, target] {
            MPCTIO tio(mpcio, thread_num);
            run_coroutines(tio, [&tio, depth, target] (yield_t &yield) {
                address_t size = address_t(1)<<depth;
                RegAS tshare;
                if (tio.player() == 2) {
                    // Send shares of the target to the computational
                    // players
                    RegAS tshare0, tshare1;
                    tshare0.randomize();
                    tshare1.set(target-tshare0.share());
                    tio.iostream_p0() << tshare0;
                    tio.iostream_p1() << tshare1;
                    printf("Using target = %016lx\n", target);
                    yield();
                } else {
                    // Get the share of the target
                    tio.iostream_server() >> tshare;
                }

                // Create a random database and sort it
                // size_t &aes_ops = tio.aes_ops();
                Duoram<RegAS> oram(tio.player(), size);
                auto A = oram.flat(tio, yield);
                A.explicitonly(true);
                // Initialize the memory to random values in parallel
                std::vector<coro_t> coroutines;
                for (address_t i=0; i<size; ++i) {
                    coroutines.emplace_back(
                        [&A, i](yield_t &yield) {
                            auto Acoro = A.context(yield);
                            RegAS v;
                            v.randomize(62);
                            Acoro[i] += v;
                        });
                }
                run_coroutines(yield, coroutines);
                A.bitonic_sort(0, depth);

                // Binary search for the target
                RegAS tindex = A.obliv_binary_search(tshare);

                // Check the answer
                if (tio.player() == 1) {
                    tio.iostream_peer() << tindex;
                } else if (tio.player() == 0) {
                    RegAS peer_tindex;
                    tio.iostream_peer() >> peer_tindex;
                    tindex += peer_tindex;
                }
                if (depth <= 10) {
                    auto check = A.reconstruct();
                    if (tio.player() == 0) {
                        for (address_t i=0;i<size;++i) {
                            printf("%04x %016lx\n", i, check[i].share());
                        }
                    }
                }
                if (tio.player() == 0) {
                    printf("Found index = %lx\n", tindex.share());
                }
            });
        });
    }
    pool.join();
}

void online_main(MPCIO &mpcio, const PRACOptions &opts, char **args)
{
    MPCTIO tio(mpcio, 0);
    if (!*args) {
        std::cerr << "Mode is required as the first argument when not preprocessing.\n";
        return;
    } else if (!strcmp(*args, "test")) {
        ++args;
        online_test(mpcio, opts, args);
    } else if (!strcmp(*args, "lamporttest")) {
        ++args;
        lamport_test(mpcio, opts, args);
    } else if (!strcmp(*args, "rdpftest")) {
        ++args;
        rdpf_test(mpcio, opts, args);
    } else if (!strcmp(*args, "rdpftime")) {
        ++args;
        rdpf_timing(mpcio, opts, args);
    } else if (!strcmp(*args, "evaltime")) {
        ++args;
        rdpfeval_timing(mpcio, opts, args);
    } else if (!strcmp(*args, "parevaltime")) {
        ++args;
        par_rdpfeval_timing(mpcio, opts, args);
    } else if (!strcmp(*args, "tupletime")) {
        ++args;
        tupleeval_timing(mpcio, opts, args);
    } else if (!strcmp(*args, "partupletime")) {
        ++args;
        par_tupleeval_timing(mpcio, opts, args);
    } else if (!strcmp(*args, "duotest")) {
        ++args;
        if (opts.use_xor_db) {
            duoram_test<RegXS>(mpcio, opts, args);
        } else {
            duoram_test<RegAS>(mpcio, opts, args);
        }
    } else if (!strcmp(*args, "cdpftest")) {
        ++args;
        cdpf_test(mpcio, opts, args);
    } else if (!strcmp(*args, "cmptest")) {
        ++args;
        compare_test(mpcio, opts, args);
    } else if (!strcmp(*args, "sorttest")) {
        ++args;
        sort_test(mpcio, opts, args);
    } else if (!strcmp(*args, "bsearch")) {
        ++args;
        bsearch_test(mpcio, opts, args);
    } else {
        std::cerr << "Unknown mode " << *args << "\n";
    }
}
