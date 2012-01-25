#include "fix16_base.h"

#ifndef FIXMATH_NO_64BIT
fix16_t fix16_mul(fix16_t a, fix16_t b)
{
    int64_t product = (int64_t)a * b;
    
    #ifndef FIXMATH_NO_OVERFLOW
    // The upper 17 bits should all be the same (the sign).
    uint32_t upper = (product >> 47);
    #endif
    
    if (product < 0)
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (~upper)
            return fix16_overflow;
        #endif
        
        #ifndef FIXMATH_NO_ROUNDING
        // This adjustment is required in order to round -1/2 correctly
        product--;
        #endif
    }
    else
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (upper)
            return fix16_overflow;
        #endif
    }
    
    #ifdef FIXMATH_NO_ROUNDING
    return product >> 16;
    #else
    fix16_t result = product >> 16;
    result += (product & 0x8000) >> 15;
    
    return result;
    #endif
}

#ifdef __GNUC__
// Count leading zeros, using processor-specific instruction if available.
#define clz(x) __builtin_clz(x)
#else
static uint8_t clz(uint32_t x)
{
    uint8_t result = 0;
    while (!(x & 0xF0000000)) { result += 4; x <<= 4; }
    while (!(x & 0x80000000)) { result += 1; x <<= 1; }
    return result;
}
#endif

fix16_t fix16_div(fix16_t a, fix16_t b)
{
    // This uses a hardware 32/32 bit division multiple times, until we have
    // computed all the bits in (a<<16)/b. Usually this takes 1-3 iterations.
    
    if (b == 0)
        return fix16_min;
    
    uint32_t remainder = (a >= 0) ? a : (-a);
    uint32_t divider = (b >= 0) ? b : (-b);
    uint32_t quotient = 0;
    int bit_pos = 17;
    
    // Kick-start the division a bit.
    // This improves speed in the worst-case scenarios where N and D are large
    // It gets a lower estimate for the result by N/(D >> 17 + 1).
    if (divider & 0xFFF00000)
    {
        uint32_t shifted_div = ((divider >> 17) + 1);
        quotient = remainder / shifted_div;
        remainder -= ((uint64_t)quotient * divider) >> 17;
    }
    
    // If the divider is divisible by 2^n, take advantage of it.
    while (!(divider & 0xF) && bit_pos >= 4)
    {
        divider >>= 4;
        bit_pos -= 4;
    }
    
    while (bit_pos >= 0)
    {
        // Shift remainder as much as we can without overflowing
        int shift = clz(remainder);
        if (shift > bit_pos) shift = bit_pos;
        remainder <<= shift;
        bit_pos -= shift;
        
        uint32_t div = remainder / divider;
        remainder = remainder % divider;
        quotient += div << bit_pos;
        
        #ifndef FIXMATH_NO_OVERFLOW
        if (div & ~(0xFFFFFFFF >> bit_pos))
            return fix16_overflow;
        #endif
        
        remainder <<= 1;
        bit_pos--;
        
        if (remainder == 0)
            break;
    }
    
    #ifndef FIXMATH_NO_ROUNDING
    // Quotient is always positive so rounding is easy
    quotient++;
    #endif
    
    fix16_t result = quotient >> 1;
    
    // Figure out the sign of the result
    if ((a ^ b) & 0x80000000)
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (result == fix16_min)
            return fix16_overflow;
        #endif
        
        result = -result;
    }
    
    return result;
}

#else /* 32-bit implementations for compilers without int64_t and for 8-bit platforms. */

fix16_t fix16_mul(fix16_t a, fix16_t b)
{
    uint32_t _a = (a >= 0) ? a : (-a);
    uint32_t _b = (b >= 0) ? b : (-b);
    
    uint8_t va[4] = {_a, (_a >> 8), (_a >> 16), (_a >> 24)};
    uint8_t vb[4] = {_b, (_b >> 8), (_b >> 16), (_b >> 24)};
    
    uint32_t low = 0;
    uint32_t mid = 0;
    
    // Result column i depends on va[0..i] and vb[i..0]

    #ifndef FIXMATH_NO_OVERFLOW
    // i = 6
    if (va[3] && vb[3]) return fix16_overflow;
    #endif
    
    // i = 5
    if (va[2] && vb[3]) mid += (uint16_t)va[2] * vb[3];
    if (va[3] && vb[2]) mid += (uint16_t)va[3] * vb[2];
    mid <<= 8;
    
    // i = 4
    if (va[1] && vb[3]) mid += (uint16_t)va[1] * vb[3];
    if (va[2] && vb[2]) mid += (uint16_t)va[2] * vb[2];
    if (va[3] && vb[1]) mid += (uint16_t)va[3] * vb[1];
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (mid & 0xFF000000) return fix16_overflow;
    #endif
    mid <<= 8;
    
    // i = 3
    if (va[0] && vb[3]) mid += (uint16_t)va[0] * vb[3];
    if (va[1] && vb[2]) mid += (uint16_t)va[1] * vb[2];
    if (va[2] && vb[1]) mid += (uint16_t)va[2] * vb[1];
    if (va[3] && vb[0]) mid += (uint16_t)va[3] * vb[0];
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (mid & 0xFF000000) return fix16_overflow;
    #endif
    mid <<= 8;
    
    // i = 2
    if (va[0] && vb[2]) mid += (uint16_t)va[0] * vb[2];
    if (va[1] && vb[1]) mid += (uint16_t)va[1] * vb[1];
    if (va[2] && vb[0]) mid += (uint16_t)va[2] * vb[0];    
    
    // i = 1
    if (va[0] && vb[1]) low += (uint16_t)va[0] * vb[1];
    if (va[1] && vb[0]) low += (uint16_t)va[1] * vb[0];
    low <<= 8;
    
    // i = 0
    if (va[0] && vb[0]) low += (uint16_t)va[0] * vb[0];
    
    #ifndef FIXMATH_NO_ROUNDING
    low += 0x8000;
    #endif
    mid += (low >> 16);
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (mid & 0x80000000)
        return fix16_overflow;
    #endif
    
    fix16_t result = mid;
    
    /* Figure out the sign of result */
    if ((a ^ b) & 0x80000000)
    {
        result = -result;
    }
    
    return result;
}

fix16_t fix16_div(fix16_t a, fix16_t b)
{
    // This uses the basic binary restoring division algorithm.
    // It appears to be faster to do the whole division manually than
    // trying to compose a 64-bit divide out of 32-bit divisions on
    // platforms without hardware divide.
    
    if (b == 0)
        return fix16_min;
    
    uint32_t remainder = (a >= 0) ? a : (-a);
    uint32_t divider = (b >= 0) ? b : (-b);

    uint32_t quotient = 0;
    uint32_t bit = 0x10000;
    
    /* The algorithm requires D >= R */
    while (divider < remainder)
    {
        divider <<= 1;
        bit <<= 1;
    }
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (!bit)
        return fix16_overflow;
    #endif
    
    if (divider & 0x80000000)
    {
        // Perform one step manually to avoid overflows later.
        // We know that divider's bottom bit is 0 here.
        if (remainder >= divider)
        {
            quotient |= bit;
            remainder -= divider;
        }
        divider >>= 1;
        bit >>= 1;
    }
    
    /* Main division loop */
    while (bit && remainder)
    {
        if (remainder >= divider)
        {
            quotient |= bit;
            remainder -= divider;
        }
        
        remainder <<= 1;
        bit >>= 1;
    }   
        
    #ifndef FIXMATH_NO_ROUNDING
    if (remainder >= divider)
    {
        quotient++;
    }
    #endif
    
    fix16_t result = quotient;
    
    /* Figure out the sign of result */
    if ((a ^ b) & 0x80000000)
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (result == fix16_min)
            return fix16_overflow;
        #endif
        
        result = -result;
    }
    
    return result;
}
#endif

fix16_t fix16_sqrt(fix16_t a)
{
    // Algorithm is quite directly from
    // http://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Binary_numeral_system_.28base_2.29
    uint32_t num = a;
    uint32_t result = 0;
    uint32_t bit;
    uint8_t n;
    
    // Many numbers will be less than 15, so
    // this gives a good balance between time spent
    // in if vs. time spent in the while loop
    // when searching for the starting value.
    if (num & 0xFFF00000)
        bit = (uint32_t)1 << 30;
    else
        bit = (uint32_t)1 << 18;
    
    while (bit > num) bit >>= 2;
    
    // The main part is executed twice, in order to avoid
    // using 64 bit values in computations.
    for (n = 0; n < 2; n++)
    {
        // First we get the top 24 bits of the answer.
        while (bit)
        {
            if (num >= result + bit)
            {
                num -= result + bit;
                result = (result >> 1) + bit;
            }
            else
            {
                result = (result >> 1);
            }
            bit >>= 2;
        }
        
        if (n == 0)
        {
            // Then process it again to get the lowest 8 bits.
            if (num > 65535)
            {
                // The remainder 'num' is too large to be shifted left
                // by 16, so we have to add 1 to result manually and
                // adjust 'num' accordingly.
                // num = a - (result + 0.5)^2
                //     = num + result^2 - (result + 0.5)^2
                //     = num - result - 0.5
                num -= result;
                num = (num << 16) - (1 << 15);
                result = (result << 16) + (1 << 15);
            }
            else
            {
                num <<= 16;
                result <<= 16;
            }
            
            bit = 1 << 14;
        }
    }

#ifndef FIXMATH_NO_ROUNDING
    // Finally, if next bit would have been 1, round the result upwards.
    if (num > result)
    {
        result++;
    }
#endif
    
    return result;
}




