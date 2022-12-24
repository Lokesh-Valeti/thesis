/* Adapted from preprocessing/bitutils.h from
 * https://git-crysp.uwaterloo.ca/avadapal/duoram by Adithya Vadapalli,
 * itself adapted from code by Ryan Henry */

#ifndef __BITUTILS_HPP__
#define __BITUTILS_HPP__

#include <cstdint>
#include <x86intrin.h>  // SSE and AVX intrinsics

static const __m128i bool128_mask[2] = {
    _mm_set_epi64x(0,1),                                        // 0b00...0001
    _mm_set_epi64x(1,0)                                         // 0b00...0001 << 64
};

static const __m128i lsb128_mask[4] = {
    _mm_setzero_si128(),                                        // 0b00...0000
    _mm_set_epi64x(0,1),                                        // 0b00...0001
    _mm_set_epi64x(0,2),                                        // 0b00...0010
    _mm_set_epi64x(0,3)                                         // 0b00...0011
};

static const __m128i lsb128_mask_inv[4] = {
    _mm_set1_epi8(-1),                                          // 0b11...1111
    _mm_set_epi64x(-1,-2),                                      // 0b11...1110
    _mm_set_epi64x(-1,-3),                                      // 0b11...1101
    _mm_set_epi64x(-1,-4)                                       // 0b11...1100
};

static const __m128i if128_mask[2] = {
    _mm_setzero_si128(),                                        // 0b00...0000
    _mm_set1_epi8(-1)                                           // 0b11...1111
};

inline __m128i xor_if(const __m128i & block1, const __m128i & block2, __m128i flag)
{
    return _mm_xor_si128(block1, _mm_and_si128(block2, flag));
}

inline __m128i xor_if(const __m128i & block1, const __m128i & block2, bool flag)
{
    return _mm_xor_si128(block1, _mm_and_si128(block2, if128_mask[flag ? 1 : 0]));
}

inline uint8_t get_lsb(const __m128i & block, uint8_t bits = 0b01)
{
    __m128i vcmp = _mm_xor_si128(_mm_and_si128(block, lsb128_mask[bits]), lsb128_mask[bits]);
    return static_cast<uint8_t>(_mm_testz_si128(vcmp, vcmp));
}

inline __m128i clear_lsb(const __m128i & block, uint8_t bits = 0b01)
{
    return _mm_and_si128(block, lsb128_mask_inv[bits]);
}

inline __m128i set_lsb(const __m128i & block, const bool val = true)
{
    return _mm_or_si128(clear_lsb(block, 0b01), lsb128_mask[val ? 0b01 : 0b00]);
}

inline uint8_t parity(const __m128i & block)
{
    // TODO...
    return 0;
}

#endif
