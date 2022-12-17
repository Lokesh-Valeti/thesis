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
static inline void prg(__m128i &out, __m128i seed, bool whichchild)
{
    AES_ECB_encrypt(out, set_lsb(seed, whichchild), prgkey.k);
}

// Compute both children of node seed
static inline void prgboth(__m128i &left, __m128i &right, __m128i seed)
{
    AES_ECB_encrypt(left, set_lsb(seed, 0), prgkey.k);
    AES_ECB_encrypt(right, set_lsb(seed, 1), prgkey.k);
}

#endif
