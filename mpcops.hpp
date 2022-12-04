#ifndef __MPCOPS_HPP__
#define __MPCOPS_HPP__

#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"

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

#endif
