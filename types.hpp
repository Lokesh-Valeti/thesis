#ifndef __OBLIVDS_TYPES_HPP__
#define __OBLIVDS_TYPES_HPP__

#include <tuple>
#include <cstdint>

#ifndef VALUE_BITS
#define VALUE_BITS 64
#endif

// Values in memory are of this type

#if VALUE_BITS == 64
typedef uint64_t value_t;
#elif VALUE_BITS == 32
typedef uint32_t value_t;
#else
#error "Unsupported value of VALUE_BITS"
#endif

// A multiplication triple is a triple (X0,Y0,Z0) held by P0 (and
// correspondingly (X1,Y1,Z1) held by P1), with all values random,
// but subject to the relation that X0*Y1 + Y0*X1 = Z0+Z1

typedef std::tuple<value_t, value_t, value_t> MultTriple;

// A half-triple is (X0,Z0) held by P0 (and correspondingly (Y1,Z1) held
// by P1), with all values random, but subject to the relation that
// X0*Y1 = Z0+Z1

typedef std::tuple<value_t, value_t> HalfTriple;

#endif
