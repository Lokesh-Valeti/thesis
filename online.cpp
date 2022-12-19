#include <bsd/stdlib.h> // arc4random_buf

#include "online.hpp"
#include "mpcops.hpp"
#include "rdpf.hpp"


static void online_test(MPCIO &mpcio, int num_threads, char **args)
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

static void lamport_test(MPCIO &mpcio, int num_threads, char **args)
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

static void rdpf_test(MPCIO &mpcio, int num_threads, char **args)
{
    nbits_t depth=6;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }

    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num, depth] {
            MPCTIO tio(mpcio, thread_num);
            if (mpcio.player == 2) {
                RDPFPair dp = tio.rdpfpair(depth);
                printf("usi0 = %016lx\n", dp.dpf[0].unit_sum_inverse);
                printf("ss0  = %016lx\n", dp.dpf[0].scaled_sum.ashare);
                printf("usi1 = %016lx\n", dp.dpf[1].unit_sum_inverse);
                printf("ss1  = %016lx\n", dp.dpf[1].scaled_sum.ashare);
            } else {
                RDPFTriple dt = tio.rdpftriple(depth);
                printf("usi0 = %016lx\n", dt.dpf[0].unit_sum_inverse);
                printf("ss0  = %016lx\n", dt.dpf[0].scaled_sum.ashare);
                printf("usi1 = %016lx\n", dt.dpf[1].unit_sum_inverse);
                printf("ss1  = %016lx\n", dt.dpf[1].scaled_sum.ashare);
                printf("usi2 = %016lx\n", dt.dpf[2].unit_sum_inverse);
                printf("ss2  = %016lx\n", dt.dpf[2].scaled_sum.ashare);
            }
        });
    }
    pool.join();
}

void online_main(MPCIO &mpcio, int num_threads, char **args)
{
    if (!*args) {
        std::cerr << "Mode is required as the first argument when not preprocessing.\n";
        return;
    } else if (!strcmp(*args, "test")) {
        ++args;
        online_test(mpcio, num_threads, args);
    } else if (!strcmp(*args, "lamporttest")) {
        ++args;
        lamport_test(mpcio, num_threads, args);
    } else if (!strcmp(*args, "rdpftest")) {
        ++args;
        rdpf_test(mpcio, num_threads, args);
    } else {
        std::cerr << "Unknown mode " << *args << "\n";
    }
}
