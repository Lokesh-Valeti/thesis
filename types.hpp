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

    RegAS() : ashare(0) {}

    inline value_t share() const { return ashare; }
    inline void set(value_t s) { ashare = s; }

    // Set each side's share to a random value nbits bits long
    inline void randomize(size_t nbits = VALUE_BITS) {
        value_t mask = MASKBITS(nbits);
        arc4random_buf(&ashare, sizeof(ashare));
        ashare &= mask;
    }

    inline RegAS &operator+=(const RegAS &rhs) {
        this->ashare += rhs.ashare;
        return *this;
    }

    inline RegAS operator+(const RegAS &rhs) const {
        RegAS res = *this;
        res += rhs;
        return res;
    }

    inline RegAS &operator-=(const RegAS &rhs) {
        this->ashare -= rhs.ashare;
        return *this;
    }

    inline RegAS operator-(const RegAS &rhs) const {
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

inline value_t combine(const RegAS &A, const RegAS &B,
        nbits_t nbits = VALUE_BITS) {
    value_t mask = ~0;
    if (nbits < VALUE_BITS) {
        mask = (value_t(1)<<nbits)-1;
    }
    return (A.ashare + B.ashare) & mask;
}

// The type of a register holding a bit share
struct RegBS {
    bit_t bshare;

    RegBS() : bshare(0) {}

    inline bit_t share() const { return bshare; }
    inline void set(bit_t s) { bshare = s; }

    // Set each side's share to a random bit
    inline void randomize() {
        unsigned char randb;
        arc4random_buf(&randb, sizeof(randb));
        bshare = randb & 1;
    }

    inline RegBS &operator^=(const RegBS &rhs) {
        this->bshare ^= rhs.bshare;
        return *this;
    }

    inline RegBS operator^(const RegBS &rhs) const {
        RegBS res = *this;
        res ^= rhs;
        return res;
    }
};

// The type of a register holding an XOR share of a value
struct RegXS {
    value_t xshare;

    RegXS() : xshare(0) {}
    RegXS(const RegBS &b) { xshare = b.bshare ? ~0 : 0; }

    inline value_t share() const { return xshare; }
    inline void set(value_t s) { xshare = s; }

    // Set each side's share to a random value nbits bits long
    inline void randomize(size_t nbits = VALUE_BITS) {
        value_t mask = MASKBITS(nbits);
        arc4random_buf(&xshare, sizeof(xshare));
        xshare &= mask;
    }

    // For RegXS, + and * should be interpreted bitwise; that is, + is
    // really XOR and * is really AND.  - is also XOR (the same as +).

    // We also include actual XOR operators for convenience

    inline RegXS &operator+=(const RegXS &rhs) {
        this->xshare ^= rhs.xshare;
        return *this;
    }

    inline RegXS operator+(const RegXS &rhs) const {
        RegXS res = *this;
        res += rhs;
        return res;
    }

    inline RegXS &operator-=(const RegXS &rhs) {
        this->xshare ^= rhs.xshare;
        return *this;
    }

    inline RegXS operator-(const RegXS &rhs) const {
        RegXS res = *this;
        res -= rhs;
        return res;
    }

    inline RegXS &operator*=(value_t rhs) {
        this->xshare &= rhs;
        return *this;
    }

    inline RegXS operator*(value_t rhs) const {
        RegXS res = *this;
        res *= rhs;
        return res;
    }

    inline RegXS &operator^=(const RegXS &rhs) {
        this->xshare ^= rhs.xshare;
        return *this;
    }

    inline RegXS operator^(const RegXS &rhs) const {
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

    // Extract a bit share of bit bitnum of the XOR-shared register
    inline RegBS bit(nbits_t bitnum) const {
        RegBS bs;
        bs.bshare = !!(xshare & (value_t(1)<<bitnum));
        return bs;
    }
};

inline value_t combine(const RegXS &A, const RegXS &B,
        nbits_t nbits = VALUE_BITS) {
    value_t mask = ~0;
    if (nbits < VALUE_BITS) {
        mask = (value_t(1)<<nbits)-1;
    }
    return (A.xshare ^ B.xshare) & mask;
}

// Some useful operations on tuples of the above types

template <typename T>
std::tuple<T,T> operator+=(std::tuple<T,T> &A,
    const std::tuple<T,T> &B)
{
    std::get<0>(A) += std::get<0>(B);
    std::get<1>(A) += std::get<1>(B);
    return A;
}

template <typename T>
std::tuple<T,T> operator+=(const std::tuple<T&,T&> &A,
    const std::tuple<T,T> &B)
{
    std::get<0>(A) += std::get<0>(B);
    std::get<1>(A) += std::get<1>(B);
    return A;
}

template <typename T>
std::tuple<T,T> operator+(const std::tuple<T,T> &A,
    const std::tuple<T,T> &B)
{
    auto res = A;
    res += B;
    return res;
}

template <typename T>
std::tuple<T,T> operator-=(const std::tuple<T&,T&> &A,
    const std::tuple<T,T> &B)
{
    std::get<0>(A) -= std::get<0>(B);
    std::get<1>(A) -= std::get<1>(B);
    return A;
}

template <typename T>
std::tuple<T,T> operator-=(std::tuple<T,T> &A,
    const std::tuple<T,T> &B)
{
    std::get<0>(A) -= std::get<0>(B);
    std::get<1>(A) -= std::get<1>(B);
    return A;
}

template <typename T>
std::tuple<T,T> operator-(const std::tuple<T,T> &A,
    const std::tuple<T,T> &B)
{
    auto res = A;
    res -= B;
    return res;
}

template <typename T>
std::tuple<T,T> operator*=(const std::tuple<T&,T&> &A,
    const std::tuple<value_t,value_t> &B)
{
    std::get<0>(A) *= std::get<0>(B);
    std::get<1>(A) *= std::get<1>(B);
    return A;
}

template <typename T>
std::tuple<T,T> operator*=(std::tuple<T,T> &A,
    const std::tuple<value_t,value_t> &B)
{
    std::get<0>(A) *= std::get<0>(B);
    std::get<1>(A) *= std::get<1>(B);
    return A;
}

template <typename T>
std::tuple<T,T> operator*(const std::tuple<T,T> &A,
    const std::tuple<value_t,value_t> &B)
{
    auto res = A;
    res *= B;
    return res;
}

template <typename T>
inline std::tuple<value_t,value_t> combine(
        const std::tuple<T,T> &A, const std::tuple<T,T> &B,
        nbits_t nbits = VALUE_BITS) {
    return std::make_tuple(
        combine(std::get<0>(A), std::get<0>(B), nbits),
        combine(std::get<1>(A), std::get<1>(B), nbits));
}

template <typename T>
std::tuple<T,T,T> operator+=(const std::tuple<T&,T&,T&> &A,
    const std::tuple<T,T,T> &B)
{
    std::get<0>(A) += std::get<0>(B);
    std::get<1>(A) += std::get<1>(B);
    std::get<2>(A) += std::get<2>(B);
    return A;
}

template <typename T>
std::tuple<T,T,T> operator+=(std::tuple<T,T,T> &A,
    const std::tuple<T,T,T> &B)
{
    std::get<0>(A) += std::get<0>(B);
    std::get<1>(A) += std::get<1>(B);
    std::get<2>(A) += std::get<2>(B);
    return A;
}

template <typename T>
std::tuple<T,T,T> operator+(const std::tuple<T,T,T> &A,
    const std::tuple<T,T,T> &B)
{
    auto res = A;
    res += B;
    return res;
}

template <typename T>
std::tuple<T,T,T> operator-=(const std::tuple<T&,T&,T&> &A,
    const std::tuple<T,T,T> &B)
{
    std::get<0>(A) -= std::get<0>(B);
    std::get<1>(A) -= std::get<1>(B);
    std::get<2>(A) -= std::get<2>(B);
    return A;
}

template <typename T>
std::tuple<T,T,T> operator-=(std::tuple<T,T,T> &A,
    const std::tuple<T,T,T> &B)
{
    std::get<0>(A) -= std::get<0>(B);
    std::get<1>(A) -= std::get<1>(B);
    std::get<2>(A) -= std::get<2>(B);
    return A;
}

template <typename T>
std::tuple<T,T,T> operator-(const std::tuple<T,T,T> &A,
    const std::tuple<T,T,T> &B)
{
    auto res = A;
    res -= B;
    return res;
}

template <typename T>
std::tuple<T,T,T> operator*=(const std::tuple<T&,T&,T&> &A,
    const std::tuple<value_t,value_t,value_t> &B)
{
    std::get<0>(A) *= std::get<0>(B);
    std::get<1>(A) *= std::get<1>(B);
    std::get<2>(A) *= std::get<2>(B);
    return A;
}

template <typename T>
std::tuple<T,T,T> operator*=(std::tuple<T,T,T> &A,
    const std::tuple<value_t,value_t,value_t> &B)
{
    std::get<0>(A) *= std::get<0>(B);
    std::get<1>(A) *= std::get<1>(B);
    std::get<2>(A) *= std::get<2>(B);
    return A;
}

template <typename T>
std::tuple<T,T,T> operator*(const std::tuple<T,T,T> &A,
    const std::tuple<value_t,value_t,value_t> &B)
{
    auto res = A;
    res *= B;
    return res;
}

template <typename T>
inline std::tuple<value_t,value_t,value_t> combine(
        const std::tuple<T,T,T> &A, const std::tuple<T,T,T> &B,
        nbits_t nbits = VALUE_BITS) {
    return std::make_tuple(
        combine(std::get<0>(A), std::get<0>(B), nbits),
        combine(std::get<1>(A), std::get<1>(B), nbits),
        combine(std::get<2>(A), std::get<2>(B), nbits));
}


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
// The *Name structs are a way to get strings representing the names of
// the types as would be given to preprocessing to create them in
// advance.
struct MultTripleName { static constexpr const char *name = "t"; };

// A half-triple is (X0,Z0) held by P0 (and correspondingly (Y1,Z1) held
// by P1), with all values random, but subject to the relation that
// X0*Y1 = Z0+Z1

using HalfTriple = std::tuple<value_t, value_t>;
struct HalfTripleName { static constexpr const char *name = "h"; };

// The type of nodes in a DPF.  This must be at least as many bits as
// the security parameter, and at least twice as many bits as value_t.

using DPFnode = __m128i;

// A Select triple is a triple of (X0,Y0,Z0) where X0 is a bit and Y0
// and Z0 are DPFnodes held by P0 (and correspondingly (X1,Y1,Z1) held
// by P1), with all values random, but subject to the relation that
// (X0*Y1) ^ (Y0*X1) = Z0^Z1.  These are only used while creating RDPFs
// in the preprocessing phase, so we never need to store them.  This is
// a struct instead of a tuple for alignment reasons.

struct SelectTriple {
    bit_t X;
    DPFnode Y, Z;
};

// These are defined in rdpf.hpp, but declared here to avoid cyclic
// header dependencies.
struct RDPFPair;
struct RDPFPairName { static constexpr const char *name = "r"; };
struct RDPFTriple;
struct RDPFTripleName { static constexpr const char *name = "r"; };
struct CDPF;
struct CDPFName { static constexpr const char *name = "c"; };

// We want the I/O (using << and >>) for many classes
// to just be a common thing: write out the bytes
// straight from memory

#define DEFAULT_IO(CLASSNAME)                  \
    template <typename T>                      \
    T& operator>>(T& is, CLASSNAME &x)         \
    {                                          \
        is.read((char *)&x, sizeof(x));        \
        return is;                             \
    }                                          \
                                               \
    template <typename T>                      \
    T& operator<<(T& os, const CLASSNAME &x)   \
    {                                          \
        os.write((const char *)&x, sizeof(x)); \
        return os;                             \
    }

// Default I/O for various types

DEFAULT_IO(DPFnode)
DEFAULT_IO(RegBS)
DEFAULT_IO(RegAS)
DEFAULT_IO(RegXS)
DEFAULT_IO(MultTriple)
DEFAULT_IO(HalfTriple)

// And for pairs and triples

#define DEFAULT_TUPLE_IO(CLASSNAME)                                  \
    template <typename T>                                            \
    T& operator>>(T& is, std::tuple<CLASSNAME, CLASSNAME> &x)        \
    {                                                                \
        is >> std::get<0>(x) >> std::get<1>(x);                      \
        return is;                                                   \
    }                                                                \
                                                                     \
    template <typename T>                                            \
    T& operator<<(T& os, const std::tuple<CLASSNAME, CLASSNAME> &x)  \
    {                                                                \
        os << std::get<0>(x) << std::get<1>(x);                      \
        return os;                                                   \
    }                                                                \
                                                                     \
    template <typename T>                                            \
    T& operator>>(T& is, std::tuple<CLASSNAME, CLASSNAME, CLASSNAME> &x) \
    {                                                                \
        is >> std::get<0>(x) >> std::get<1>(x) >> std::get<2>(x);    \
        return is;                                                   \
    }                                                                \
                                                                     \
    template <typename T>                                            \
    T& operator<<(T& os, const std::tuple<CLASSNAME, CLASSNAME, CLASSNAME> &x) \
    {                                                                \
        os << std::get<0>(x) << std::get<1>(x) << std::get<2>(x);    \
        return os;                                                   \
    }

DEFAULT_TUPLE_IO(RegAS)
DEFAULT_TUPLE_IO(RegXS)

#endif
