#ifndef __OBLIVDS_TYPES_HPP__
#define __OBLIVDS_TYPES_HPP__

#include <tuple>
#include <cstdint>
#include <x86intrin.h>  // SSE and AVX intrinsics
#include <bsd/stdlib.h> // arc4random_buf

// The number of bits in an MPC secret-shared memory word

#ifndef VALUE_BITS
#define VALUE_BITS 64
#endif

// Values in MPC secret-shared memory are of this type.
// This is the type of the underlying shared value, not the types of the
// shares themselves.

#if VALUE_BITS == 64
using value_t = uint64_t;
#elif VALUE_BITS == 32
using value_t = uint32_t;
#else
#error "Unsupported value of VALUE_BITS"
#endif

// Secret-shared bits are of this type.  Note that it is standards
// compliant to treat a bool as an unsigned integer type with values 0
// and 1.

using bit_t = bool;

// Counts of the number of bits in a value are of this type, which must
// be large enough to store the _value_ VALUE_BITS
using nbits_t = uint8_t;

// Convert a number of bits to the number of bytes required to store (or
// more to the point, send) them.
#define BITBYTES(nbits) (((nbits)+7)>>3)

// A mask of this many bits; the test is to prevent 1<<nbits from
// overflowing if nbits == VALUE_BITS
#define MASKBITS(nbits) (((nbits) < VALUE_BITS) ? (value_t(1)<<(nbits))-1 : ~0)

// The type of a register holding an additive share of a value
struct RegAS {
    value_t ashare;

    // Set each side's share to a random value nbits bits long
    inline void randomize(size_t nbits = VALUE_BITS) {
        value_t mask = MASKBITS(nbits);
        arc4random_buf(&ashare, sizeof(ashare));
        ashare &= mask;
    }

    inline RegAS &operator+=(RegAS &rhs) {
        this->ashare += rhs.ashare;
        return *this;
    }

    inline RegAS operator+(RegAS &rhs) const {
        RegAS res = *this;
        res += rhs;
        return res;
    }

    inline RegAS &operator-=(RegAS &rhs) {
        this->ashare -= rhs.ashare;
        return *this;
    }

    inline RegAS operator-(RegAS &rhs) const {
        RegAS res = *this;
        res -= rhs;
        return res;
    }

    inline RegAS &operator*=(value_t rhs) {
        this->ashare *= rhs;
        return *this;
    }

    inline RegAS operator*(value_t rhs) const {
        RegAS res = *this;
        res *= rhs;
        return res;
    }

    inline RegAS &operator&=(value_t mask) {
        this->ashare &= mask;
        return *this;
    }

    inline RegAS operator&(value_t mask) const {
        RegAS res = *this;
        res &= mask;
        return res;
    }
};

// The type of a register holding an XOR share of a value
struct RegXS {
    value_t xshare;

    // Set each side's share to a random value nbits bits long
    inline void randomize(size_t nbits = VALUE_BITS) {
        value_t mask = MASKBITS(nbits);
        arc4random_buf(&xshare, sizeof(xshare));
        xshare &= mask;
    }

    inline RegXS &operator^=(RegXS &rhs) {
        this->xshare ^= rhs.xshare;
        return *this;
    }

    inline RegXS operator^(RegXS &rhs) const {
        RegXS res = *this;
        res ^= rhs;
        return res;
    }

    inline RegXS &operator&=(value_t mask) {
        this->xshare &= mask;
        return *this;
    }

    inline RegXS operator&(value_t mask) const {
        RegXS res = *this;
        res &= mask;
        return res;
    }
};

// The type of a register holding a bit share
struct RegBS {
    bit_t bshare;

    // Set each side's share to a random bit
    inline void randomize() {
        arc4random_buf(&bshare, sizeof(bshare));
        bshare &= 1;
    }
};

// The _maximum_ number of bits in an MPC address; the actual size of
// the memory will typically be set at runtime, but it cannot exceed
// this value.  It is more efficient (in terms of communication) in some
// places for this value to be at most 32.

#ifndef ADDRESS_MAX_BITS
#define ADDRESS_MAX_BITS 32
#endif

// Addresses of MPC secret-shared memory are of this type

#if ADDRESS_MAX_BITS <= 32
using address_t = uint32_t;
#elif ADDRESS_MAX_BITS <= 64
using address_t = uint64_t;
#else
#error "Unsupported value of ADDRESS_MAX_BITS"
#endif

#if ADDRESS_MAX_BITS > VALUE_BITS
#error "VALUE_BITS must be at least as large as ADDRESS_MAX_BITS"
#endif

// A multiplication triple is a triple (X0,Y0,Z0) held by P0 (and
// correspondingly (X1,Y1,Z1) held by P1), with all values random,
// but subject to the relation that X0*Y1 + Y0*X1 = Z0+Z1

using MultTriple = std::tuple<value_t, value_t, value_t>;

// A half-triple is (X0,Z0) held by P0 (and correspondingly (Y1,Z1) held
// by P1), with all values random, but subject to the relation that
// X0*Y1 = Z0+Z1

using HalfTriple = std::tuple<value_t, value_t>;

// The type of nodes in a DPF.  This must be at least as many bits as
// the security parameter, and at least twice as many bits as value_t.

using DPFnode = __m128i;

// An AND triple is a triple of (X0,Y0,Z0) of DPFnodes held by P0 (and
// correspondingly (X1,Y1,Z1) held by P1), with all values random, but
// subject to the relation that (X0&Y1) ^ (Y0&X1) = Z0^Z1.  These are
// only used while creating RDPFs in the preprocessing phase, so we
// never need to store them.

using AndTriple = std::tuple<DPFnode, DPFnode, DPFnode>;

#endif
