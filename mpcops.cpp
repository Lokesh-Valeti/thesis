#include "mpcops.hpp"

// as_ denotes additive shares
// xs_ denotes xor shares
// bs_ denotes a share of a single bit (which is effectively both an xor
//     share and an additive share mod 2)

// P0 and P1 both hold additive shares of x (shares are x0 and x1) and y
// (shares are y0 and y1); compute additive shares of z = x*y =
// (x0+x1)*(y0+y1). x, y, and z are each at most nbits bits long.
//
// Cost:
// 2 words sent in 1 message
// consumes 1 MultTriple
void mpc_mul(MPCTIO &tio, yield_t &yield,
    value_t &as_z, value_t as_x, value_t as_y,
    nbits_t nbits)
{
    const value_t mask = MASKBITS(nbits);
    // Compute as_z to be an additive share of (x0*y1+y0*x1)
    mpc_cross(tio, yield, as_z, as_x, as_y, nbits);
    // Add x0*y0 (the peer will add x1*y1)
    as_z = (as_z + as_x * as_y) & mask;
}

// P0 and P1 both hold additive shares of x (shares are x0 and x1) and y
// (shares are y0 and y1); compute additive shares of z = x0*y1 + y0*x1.
// x, y, and z are each at most nbits bits long.
//
// Cost:
// 2 words sent in 1 message
// consumes 1 MultTriple
void mpc_cross(MPCTIO &tio, yield_t &yield,
    value_t &as_z, value_t as_x, value_t as_y,
    nbits_t nbits)
{
    const value_t mask = MASKBITS(nbits);
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

    as_z = ((as_x * peer_blind_y) - (Y * peer_blind_x) + Z) & mask;
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
    const value_t mask = MASKBITS(nbits);
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

// P0 and P1 hold bit shares f0 and f1 of the single bit f, and additive
// shares y0 and y1 of the value y; compute additive shares of
// z = f * y = (f0 XOR f1) * (y0 + y1).  y and z are each at most nbits
// bits long.
//
// Cost:
// 2 words sent in 1 message
// consumes 1 MultTriple
void mpc_flagmult(MPCTIO &tio, yield_t &yield,
    value_t &as_z, bit_t bs_f, value_t as_y,
    nbits_t nbits)
{
    const value_t mask = MASKBITS(nbits);

    // Compute additive shares of [(1-2*f0)*y0]*f1 + [(1-2*f1)*y1]*f0
    value_t bs_fval = value_t(bs_f);
    mpc_cross(tio, yield, as_z, (1-2*bs_fval)*as_y, bs_fval, nbits);

    // Add f0*y0 (and the peer will add f1*y1)
    as_z = (as_z + bs_fval*as_y) & mask;

    // Now the shares add up to:
    // [(1-2*f0)*y0]*f1 + [(1-2*f1)*y1]*f0 + f0*y0 + f1*y1
    // which you can rearrange to see that it's equal to the desired
    // (f0 + f1 - 2*f0*f1)*(y0+y1), since f0 XOR f1 = (f0 + f1 - 2*f0*f1).
}

// P0 and P1 hold bit shares f0 and f1 of the single bit f, and additive
// shares of the values x and y; compute additive shares of z, where
// z = x if f=0 and z = y if f=1.  x, y, and z are each at most nbits
// bits long.
//
// Cost:
// 2 words sent in 1 message
// consumes 1 MultTriple
void mpc_select(MPCTIO &tio, yield_t &yield,
    value_t &as_z, bit_t bs_f, value_t as_x, value_t as_y,
    nbits_t nbits)
{
    const value_t mask = MASKBITS(nbits);

    // The desired result is z = x + f * (y-x)
    mpc_flagmult(tio, yield, as_z, bs_f, as_y-as_x, nbits);
    as_z = (as_z + as_x) & mask;
}

// P0 and P1 hold bit shares f0 and f1 of the single bit f, and additive
// shares of the values x and y. Obliviously swap x and y; that is,
// replace x and y with new additive sharings of x and y respectively
// (if f=0) or y and x respectively (if f=1).  x and y are each at most
// nbits bits long.
//
// Cost:
// 2 words sent in 1 message
// consumes 1 MultTriple
void mpc_oswap(MPCTIO &tio, yield_t &yield,
    value_t &as_x, value_t &as_y, bit_t bs_f,
    nbits_t nbits)
{
    const value_t mask = MASKBITS(nbits);

    // Let s = f*(y-x).  Then the desired result is
    // x <- x + s, y <- y - s.
    value_t as_s;
    mpc_flagmult(tio, yield, as_s, bs_f, as_y-as_x, nbits);
    as_x = (as_x + as_s) & mask;
    as_y = (as_y - as_s) & mask;
}

// P0 and P1 hold XOR shares of x. Compute additive shares of the same
// x. x is at most nbits bits long.
//
// Cost:
// nbits-1 words sent in 1 message
// consumes nbits-1 HalfTriples
void mpc_xs_to_as(MPCTIO &tio, yield_t &yield,
    value_t &as_x, value_t xs_x,
    nbits_t nbits)
{
    const value_t mask = MASKBITS(nbits);

    // We use the fact that for any nbits-bit A and B,
    // A+B = (A XOR B) + 2*(A AND B)  mod 2^nbits
    // so if we have additive shares C0 and C1 of 2*(A AND B)
    // (so C0 + C1 = 2*(A AND B)), then (A-C0) and (B-C1) are
    // additive shares of (A XOR B).

    // To get additive shares of 2*(A AND B) (mod 2^nbits), we first
    // note that we can ignore the top bits of A and B, since the
    // multiplication by 2 will shift it out of the nbits-bit range.
    // For the other bits, use valuemult to get the product of the
    // corresponding bit i of A and B (i=0..nbits-2), and compute
    // C = \sum_i 2^{i+1} * (A_i * B_i).

    // This can all be done in a single message, using the coroutine
    // mechanism to have all nbits-1 instances of valuemult queue their
    // message, then yield, so that all of their messages get sent at
    // once, then each will read their results.

    value_t as_bitand[nbits-1];
    std::vector<coro_t> coroutines;
    for (nbits_t i=0; i<nbits-1; ++i) {
        coroutines.emplace_back(
            [&](yield_t &yield) {
                mpc_valuemul(tio, yield, as_bitand[i], (xs_x>>i)&1, nbits);
            });
    }
    run_coroutines(yield, coroutines);
    value_t as_C = 0;
    for (nbits_t i=0; i<nbits-1; ++i) {
        as_C += (as_bitand[i]<<(i+1));
    }
    as_x = (xs_x - as_C) & mask;
}
