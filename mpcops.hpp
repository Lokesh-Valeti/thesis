#ifndef __MPCOPS_HPP__
#define __MPCOPS_HPP__

#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"

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
    nbits_t nbits = VALUE_BITS);

// P0 and P1 both hold additive shares of x (shares are x0 and x1) and y
// (shares are y0 and y1); compute additive shares of z = x0*y1 + y0*x1.
// x, y, and z are each at most nbits bits long.
//
// Cost:
// 2 words sent in 1 message
// consumes 1 MultTriple
void mpc_cross(MPCTIO &tio, yield_t &yield,
    value_t &as_z, value_t as_x, value_t as_y,
    nbits_t nbits = VALUE_BITS);

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
    nbits_t nbits = VALUE_BITS);

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
    nbits_t nbits = VALUE_BITS);

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
    nbits_t nbits = VALUE_BITS);

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
    nbits_t nbits = VALUE_BITS);

// P0 and P1 hold XOR shares of x. Compute additive shares of the same
// x. x is at most nbits bits long.
//
// Cost:
// nbits-1 words sent in 1 message
// consumes nbits-1 HalfTriples
void mpc_xs_to_as(MPCTIO &tio, yield_t &yield,
    value_t &as_x, value_t xs_x,
    nbits_t nbits = VALUE_BITS);

#endif
