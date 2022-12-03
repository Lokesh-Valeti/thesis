#include <bsd/stdlib.h> // arc4random_buf

#include "online.hpp"
#include "mpcops.hpp"


void online_comp(MPCIO &mpcio, int num_threads, char **args)
{
    nbits_t nbits = VALUE_BITS;

    if (*args) {
        nbits = atoi(*args);
    }

    value_t A[5];

    arc4random_buf(A, 3*sizeof(value_t));
    std::cout << A[0] << "\n";
    std::cout << A[1] << "\n";
    std::cout << A[2] << "\n";
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(
        [&](yield_t &yield) {
            MultTriple T;
            mpcio.triples[0].get(T);
            mpc_mul(mpcio, 0, yield, A[3], A[0], A[1], T, nbits);
        });
    coroutines.emplace_back(
        [&](yield_t &yield) {
            HalfTriple H;
            mpcio.halftriples[0].get(H);
            mpc_valuemul(mpcio, 0, yield, A[4], A[2], H, nbits);
        });
    run_coroutines(mpcio, coroutines);
    std::cout << A[3] << "\n";
    std::cout << A[4] << "\n";
}

void online_server(MPCServerIO &mpcio, char **args)
{
}
