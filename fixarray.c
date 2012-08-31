#include "fixarray.h"
#include <string.h> /* For memcpy() */

#ifdef FIXMATH_NO_64BIT

// Calculates the dotproduct of two vectors of size n.
// If overflow happens, sets flag in errors
fix16_t fa16_dot(const fix16_t *a, uint_fast8_t a_stride,
                 const fix16_t *b, uint_fast8_t b_stride,
                 uint_fast8_t n)
{
    fix16_t sum = 0;
    
    while (n--)
    {
        // Compute result
        if (*a != 0 && *b != 0)
        {
            fix16_t product = fix16_mul(*a, *b);
            sum = fix16_add(sum, product);
            
            if (sum == fix16_overflow || product == fix16_overflow)
                return fix16_overflow;
        }
        
        // Go to next item
        a += a_stride;
        b += b_stride;
    }
    
    return sum;
}

#else

// Because dotproduct() is the hotspot of matrix multiplication,
// it has a specialized 64-bit routine in addition to the normal
// fix16_mul()-based one. This is especially efficient on ARM processors
// which have SMLAL instruction.
fix16_t fa16_dot(const fix16_t *a, uint_fast8_t a_stride,
                 const fix16_t *b, uint_fast8_t b_stride,
                 uint_fast8_t n)
{
    int64_t sum = 0;
    
    while (n--)
    {
        if (*a != 0 && *b != 0)
        {
            sum += (int64_t)(*a) * (*b);
        }
        
        // Go to next item
        a += a_stride;
        b += b_stride;
    }
    
    // The upper 17 bits should all be the same (the sign).
    uint32_t upper = sum >> 47;
    if (sum < 0)
    {
        upper = ~upper;
        
        #ifndef FIXMATH_NO_ROUNDING
        // This adjustment is required in order to round -1/2 correctly
        sum--;
        #endif
    }
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (upper)
        return fix16_overflow;
    #endif
    
    fix16_t result = sum >> 16;

    #ifndef FIXMATH_NO_ROUNDING
    result += (sum & 0x8000) >> 15;
    #endif
    
    return result;
}
#endif

#ifdef __GNUC__
// Count leading zeros, using processor-specific instruction if available.
#define clz(x) (__builtin_clzl(x) - (8 * sizeof(long) - 32))
#else
static uint8_t clz(uint32_t x)
{
  uint8_t result = 0;
  if (x == 0) return 32;
  while (!(x & 0xF0000000)) { result += 4; x <<= 4; }
  while (!(x & 0x80000000)) { result += 1; x <<= 1; }
  return result;
}
#endif

static fix16_t scale_value(fix16_t value, int_fast8_t scale)
{
    if (scale > 0)
    {
        fix16_t temp = value << scale;
        if (temp >> scale != value)
            return fix16_overflow;
        else
            return temp;
    }
    else if (scale < 0)
    {
        return value >> -scale;
    }
    else
    {
        return value;
    }
}

#ifndef FIXMATH_NO_64BIT
// Calculates the norm of a vector
fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n)
{
    int64_t sum = 0;
    
    while (n--)
    {
        if (*a != 0)
        {
            sum += (int64_t)(*a) * (*a);
        }
        
        a += a_stride;
    }
    
    int_fast8_t scale = 0;
    uint32_t highpart = (uint32_t)(sum >> 32);
    uint32_t lowpart = (uint32_t)sum;
    if (highpart)
        scale = 33 - clz(highpart);
    else if (lowpart & 0x80000000)
        scale = 1;
    
    if (scale & 1) scale++;
    
    fix16_t result = fix16_sqrt((uint32_t)(sum >> scale));
    result = scale_value(result, scale / 2 - 8);
    
    return result;
}

#else

static uint_fast8_t ilog2(uint_fast8_t v)
{
    uint_fast8_t result = 0;
    if (v & 0xF0) { result += 4; v >>= 4; }
    while (v) { result++; v >>= 1; }
    return result;
}

fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n)
{
    fix16_t sum = 0;
    fix16_t max = 0;
    
    // Calculate inclusive OR of all values to find out the maximum.
    {
        uint_fast8_t i;
        const fix16_t *p = a;
        for (i = 0; i < n; i++, p += a_stride)
        {
            max |= fix16_abs(*p);
        }
    }
    
    // To avoid overflows, the values before squaring can be max 128.0,
    // i.e. v & 0xFF800000 must be 0. Also, to avoid overflow in sum,
    // we need additional log2(n) bits of space.
    int_fast8_t scale = clz(max) - 9 - ilog2(n) / 2;
    
    while (n--)
    {
        fix16_t val = scale_value(*a, scale);
        fix16_t product = fix16_mul(val, val);
        sum = fix16_add(sum, product);
        
        a += a_stride;
    }
    
    if (sum == fix16_overflow)
        return sum;
    
    fix16_t result = fix16_sqrt(sum);
    return scale_value(result, -scale);
}

#endif

void fa16_unalias(void *dest, void **a, void **b, void *tmp, unsigned size)
{
    if (dest == *a)
    {
        memcpy(tmp, *a, size);
        *a = tmp;
        
        if (dest == *b)
            *b = tmp;
    }
    else if (dest == *b)
    {
        memcpy(tmp, *b, size);
        *b = tmp;
    }
}

