#include "mpcops.hpp"

// as_ denotes additive shares
// xs_ denotes xor shares
// bs_ denotes a share of a single bit (which is effectively both an xor
//     share and an additive share mod 2)

// P0 and P1 both hold additive shares of x and y; compute additive
// shares of z = x*y. x, y, and z are each at most nbits bits long.
//
// Cost:
// 1 word sent in 1 message
// consumes 1 MultTriple
void mpc_mul(MPCTIO &tio, yield_t &yield,
    value_t &as_z, value_t as_x, value_t as_y,
    nbits_t nbits)
{
    value_t mask = MASKBITS(nbits);
    size_t nbytes = BITBYTES(nbits);
    auto [X, Y, Z] = tio.triple();

    // Send x+X and y+Y
    value_t blind_x = (as_x + X) & mask;
    value_t blind_y = (as_y + Y) & mask;

    tio.queue_peer(&blind_x, nbytes);
    tio.queue_peer(&blind_y, nbytes);

    yield();

    // Read the peer's x+X and y+Y
    value_t  peer_blind_x, peer_blind_y;
    tio.recv_peer(&peer_blind_x, nbytes);
    tio.recv_peer(&peer_blind_y, nbytes);

    as_z = ((as_x * (as_y + peer_blind_y)) - Y * peer_blind_x + Z) & mask;
}

// P0 holds the (complete) value x, P1 holds the (complete) value y;
// compute additive shares of z = x*y.  x, y, and z are each at most
// nbits bits long.  The parameter is called x, but P1 will pass y
// there.
//
// Cost:
// 1 word sent in 1 message
// consumes 1 HalfTriple
void mpc_valuemul(MPCTIO &tio, yield_t &yield,
    value_t &as_z, value_t x,
    nbits_t nbits)
{
    value_t mask = MASKBITS(nbits);
    size_t nbytes = BITBYTES(nbits);
    auto [X, Z] = tio.halftriple();

    // Send x+X
    value_t blind_x = (x + X) & mask;

    tio.queue_peer(&blind_x, nbytes);

    yield();

    // Read the peer's y+Y
    value_t  peer_blind_y;
    tio.recv_peer(&peer_blind_y, nbytes);

    if (tio.player() == 0) {
        as_z = ((x * peer_blind_y) + Z) & mask;
    } else if (tio.player() == 1) {
        as_z = ((-X * peer_blind_y) + Z) & mask;
    }
}
