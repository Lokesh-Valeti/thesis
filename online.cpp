#include <bsd/stdlib.h> // arc4random_buf

#include "online.hpp"
#include "mpcops.hpp"


void online_comp(MPCIO &mpcio, int num_threads, char **args)
{
    nbits_t nbits = VALUE_BITS;

    if (*args) {
        nbits = atoi(*args);
    }

    size_t memsize = 13;

    MPCTIO tio(mpcio, 0);

    value_t *A = new value_t[memsize];

    arc4random_buf(A, memsize*sizeof(value_t));
    A[5] &= 1;
    A[8] &= 1;
    printf("A:\n"); for (size_t i=0; i<memsize; ++i) printf("%3lu: %016lX\n", i, A[i]);
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
    printf("\n");
    printf("A:\n"); for (size_t i=0; i<memsize; ++i) printf("%3lu: %016lX\n", i, A[i]);

    // Check the answers
    if (mpcio.player) {
        tio.queue_peer(A, memsize*sizeof(value_t));
        tio.send();
    } else {
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

void online_server(MPCServerIO &mpcio, int num_threads, char **args)
{
}
