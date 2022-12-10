#include <bsd/stdlib.h> // arc4random_buf

#include "online.hpp"
#include "mpcops.hpp"


static void online_test(MPCIO &mpcio, int num_threads, char **args)
{
    nbits_t nbits = VALUE_BITS;

    if (*args) {
        nbits = atoi(*args);
    }

    size_t memsize = 13;

    MPCTIO tio(mpcio, 0);
    bool is_server = (mpcio.player == 2);

    value_t *A = new value_t[memsize];

    if (!is_server) {
        arc4random_buf(A, memsize*sizeof(value_t));
        A[5] &= 1;
        A[8] &= 1;
        printf("A:\n"); for (size_t i=0; i<memsize; ++i) printf("%3lu: %016lX\n", i, A[i]);
    }
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_mul(tio, yield, A[2], A[0], A[1], nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_valuemul(tio, yield, A[4], A[3], nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_flagmult(tio, yield, A[7], A[5], A[6], nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_oswap(tio, yield, A[9], A[10], A[8], nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_xs_to_as(tio, yield, A[12], A[11], nbits);
        });
    run_coroutines(tio, coroutines);
    if (!is_server) {
        printf("\n");
        printf("A:\n"); for (size_t i=0; i<memsize; ++i) printf("%3lu: %016lX\n", i, A[i]);
    }

    // Check the answers
    if (mpcio.player == 1) {
        tio.queue_peer(A, memsize*sizeof(value_t));
        tio.send();
    } else if (mpcio.player == 0) {
        value_t *B = new value_t[memsize];
        value_t *S = new value_t[memsize];
        tio.recv_peer(B, memsize*sizeof(value_t));
        for(size_t i=0; i<memsize; ++i) S[i] = A[i]+B[i];
        printf("S:\n"); for (size_t i=0; i<memsize; ++i) printf("%3lu: %016lX\n", i, S[i]);
        printf("\n%016lx\n", S[0]*S[1]-S[2]);
        printf("%016lx\n", (A[3]*B[3])-S[4]);
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
    } else {
        std::cerr << "Unknown mode " << *args << "\n";
    }
}
