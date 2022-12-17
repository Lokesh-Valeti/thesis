#ifndef __PRG_HPP__
#define __PRG_HPP__

#include "bitutils.hpp"
#include "aes.hpp"

static const struct PRGkey {
    AESkey k;
    PRGkey(__m128i key = _mm_set_epi64x(314159265, 271828182)) {
        AES_128_Key_Expansion(k, key);
    }
} prgkey;

// Compute one of the children of node seed; whichchild=0 for
// the left child, 1 for the right child
static inline void prg(__m128i &out, __m128i seed, bool whichchild,
    size_t &op_counter)
{
    __m128i in = set_lsb(seed, whichchild);
    __m128i mid;
    AES_ECB_encrypt(mid, set_lsb(seed, whichchild), prgkey.k, op_counter);
    out = _mm_xor_si128(mid, in);
}

// Compute both children of node seed
static inline void prgboth(__m128i &left, __m128i &right, __m128i seed,
    size_t &op_counter)
{
    __m128i in0 = set_lsb(seed, 0);
    __m128i in1 = set_lsb(seed, 1);
    __m128i mid0, mid1;
    AES_ECB_encrypt(mid0, set_lsb(seed, 0), prgkey.k, op_counter);
    AES_ECB_encrypt(mid1, set_lsb(seed, 1), prgkey.k, op_counter);
    left = _mm_xor_si128(mid0, in0);
    right = _mm_xor_si128(mid1, in1);
}

#endif
