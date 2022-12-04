#include <bsd/stdlib.h> // arc4random_buf

#include "online.hpp"
#include "mpcops.hpp"


void online_comp(MPCIO &mpcio, int num_threads, char **args)
{
    nbits_t nbits = VALUE_BITS;

    if (*args) {
        nbits = atoi(*args);
    }

    size_t memsize = 5;

    MPCTIO tio(mpcio, 0);

    value_t *A = new value_t[memsize];

    arc4random_buf(A, 3*sizeof(value_t));
    std::cout << A[0] << "\n";
    std::cout << A[1] << "\n";
    std::cout << A[2] << "\n";
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_mul(tio, yield, A[3], A[0], A[1], nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            mpc_valuemul(tio, yield, A[4], A[2], nbits);
        });
    run_coroutines(tio, coroutines);
    std::cout << A[3] << "\n";
    std::cout << A[4] << "\n";

    // Check the answers
    if (mpcio.player) {
        tio.queue_peer(A, memsize*sizeof(value_t));
        tio.send();
    } else {
        value_t *B = new value_t[memsize];
        tio.recv_peer(B, memsize*sizeof(value_t));
        printf("%016lx\n", ((A[0]+B[0])*(A[1]+B[1])-(A[3]+B[3])));
        printf("%016lx\n", (A[2]*B[2])-(A[4]+B[4]));
        delete[] B;
    }

    delete[] A;
}

void online_server(MPCServerIO &mpcio, int num_threads, char **args)
{
}
