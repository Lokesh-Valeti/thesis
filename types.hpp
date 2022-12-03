#ifndef __OBLIVDS_TYPES_HPP__
#define __OBLIVDS_TYPES_HPP__

#include <tuple>
#include <cstdint>

// The number of bits in an MPC secret-shared memory word

#ifndef VALUE_BITS
#define VALUE_BITS 64
#endif

// Values in MPC secret-shared memory are of this type

#if VALUE_BITS == 64
typedef uint64_t value_t;
#elif VALUE_BITS == 32
typedef uint32_t value_t;
#else
#error "Unsupported value of VALUE_BITS"
#endif

// The _maximum_ number of bits in an MPC address; the actual size of
// the memory will typically be set at runtime, but it cannot exceed
// this value.  It is more efficient (in terms of communication) in some
// places for this value to be at most 32.

#ifndef ADDRESS_MAX_BITS
#define ADDRESS_MAX_BITS 32
#endif

// Addresses of MPC secret-shared memory are of this type

#if ADDRESS_MAX_BITS <= 32
typedef uint32_t address_t;
#elif ADDRESS_MAX_BITS <= 64
typedef uint64_t address_t;
#else
#error "Unsupported value of ADDRESS_MAX_BITS"
#endif

#if ADDRESS_MAX_BITS > VALUE_BITS
#error "VALUE_BITS must be at least as large as ADDRESS_MAX_BITS"
#endif

// Secret-shared bits are of this type.  Note that it is standards
// compliant to treat a bool as an unsigned integer type with values 0
// and 1.

typedef bool bit_t;

// A multiplication triple is a triple (X0,Y0,Z0) held by P0 (and
// correspondingly (X1,Y1,Z1) held by P1), with all values random,
// but subject to the relation that X0*Y1 + Y0*X1 = Z0+Z1

typedef std::tuple<value_t, value_t, value_t> MultTriple;

// A half-triple is (X0,Z0) held by P0 (and correspondingly (Y1,Z1) held
// by P1), with all values random, but subject to the relation that
// X0*Y1 = Z0+Z1

typedef std::tuple<value_t, value_t> HalfTriple;

#endif
