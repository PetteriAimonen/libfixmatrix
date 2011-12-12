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

fix16_t fix16_div(fix16_t a, fix16_t b)
{
    #ifdef FIXMATH_NO_ROUNDING
    int64_t quotient = ((int64_t)a << 16) / b;
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (quotient >> 63 != quotient >> 31)
        return fix16_overflow;
    #endif
    
    return quotient;
    #else
    // To implement rounding, we first shift temp left by 17 bits instead of
    // the 16 required by the fixed point format. The lowest bit is used
    // for rounding.
    // a/b +- 0.5 = (2a/b +- 1)/2
    int64_t a_x17 = (int64_t)a << 17;
    int64_t quotient = a_x17 / b;

    #ifndef FIXMATH_NO_OVERFLOW
    if (quotient >> 63 != quotient >> 32)
        return fix16_overflow;
    #endif
    
    // Now is the time to subtract 1 for negative quotient and
    // add 1 for positive quotient. However we will do it after
    // "dividing" by 2, because right shift is not really same
    // as dividing for negative numbers. This way we can fix
    // two things at once.
    fix16_t result = quotient >> 1;
    
    // Now depending on result sign and quotient lowest bit add:
    // negative, 0:   0 + 0 =  0
    // negative, 1:  -1 + 1 =  0
    // positive, 0:   0 + 0 =  0
    // positive, 1:   1 + 0 =  1
    //                ^   ^--- Fix for shift vs. divide
    //                `------- Fix for quotient +-= 1 
    if (result >= 0) result += quotient & 1;
    
    return result;
    #endif
}

#else /* 32-bit implementations for compilers without int64_t. */

fix16_t fix16_mul(fix16_t a, fix16_t b)
{
    // Each argument is divided to 16-bit parts.
    //          AB
    //      *   CD
    // -----------
    //          BD  16 * 16 -> 32 bit products
    //         CB
    //         AD
    //        AC
    //       |----| 64 bit product
    int32_t A = (a >> 16), C = (b >> 16);
    uint32_t B = (a & 0xFFFF), D = (b & 0xFFFF);
    
    int32_t AC = A*C;
    int32_t AD_CB = A*D + C*B;
    uint32_t BD = B*D;
    
    int32_t product_hi = AC + (AD_CB >> 16);
    
    // Handle carry from lower 32 bits to upper part of result.
    uint32_t ad_cb_temp = AD_CB << 16;
    uint32_t product_lo = BD + ad_cb_temp;
    if (product_lo < BD)
        product_hi++;
    
#ifndef FIXMATH_NO_OVERFLOW
    // The upper 17 bits should all be the same (the sign).
    if (product_hi >> 31 != product_hi >> 15)
        return fix16_overflow;
#endif
    
#ifdef FIXMATH_NO_ROUNDING
    return (product_hi << 16) | (product_lo >> 16);
#else
    // Subtracting 0x8000 (= 0.5) and then using signed right shift
    // achieves proper rounding to result-1, except in the corner
    // case of negative numbers and lowest word = 0x8000.
    // To handle that, we also have to subtract 1 for negative numbers.
    uint32_t product_lo_tmp = product_lo;
    product_lo -= 0x8000;
    product_lo -= (uint32_t)product_hi >> 31;
    if (product_lo > product_lo_tmp)
        product_hi--;
    
    // Discard the lowest 16 bits. Note that this is not exactly the same
    // as dividing by 0x10000. For example if product = -1, result will
    // also be -1 and not 0. This is compensated by adding +1 to the result
    // and compensating this in turn in the rounding above.
    fix16_t result = (product_hi << 16) | (product_lo >> 16);
    result += 1;
    return result;
#endif
}

#ifdef __GNUC__
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

#include <stdio.h>
fix16_t fix16_div(fix16_t a, fix16_t b)
{
    // This uses the basic binary restoring division algorithm.
    // It appears to be faster to do the whole division manually than
    // trying to compose a 64-bit divide out of 32-bit divisions and
    // multiplications, especially on platforms without hardware divide.
    
    if (b == 0)
        return fix16_min;
    
    uint32_t remainder = (a >= 0) ? a : (-a);
    uint32_t divider = (b >= 0) ? b : (-b);
    uint32_t quotient = 0;
    int bit_pos = 17;
    
//     printf("--\n");
//     printf("%12u/%12u %08x %d\n", remainder, divider, quotient, bit_pos);
    
    // Kick-start the division a bit.
    // This improves speed in the worst-case scenarios where N and D are large
    // It gets a lower estimate for the result by N/(D >> 17 + 1).
    if (divider & 0xFFF00000)
    {
        printf("here\n");
        uint32_t shifted_div = ((divider >> 17) + 1);
        quotient = remainder / shifted_div;
        remainder -= ((uint64_t)quotient * divider) >> 17;
    }
    
//     printf("%12u/%12u %08x %d\n", remainder, divider, quotient, bit_pos);
    
    while (!(divider & 1) && bit_pos >= 1)
    {
        divider >>= 1;
        bit_pos -= 1;
    }
    
    while (bit_pos >= 0)
    {
        int shift = clz(remainder);
        if (shift > bit_pos) shift = bit_pos;
        remainder <<= shift;
        bit_pos -= shift;
        
//         printf("%12u/%12u %08x %d\n", remainder, divider, quotient, bit_pos);
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
    
//     printf("%12u/%12u %08x %d\n", remainder, divider, quotient, bit_pos);
    
    #ifndef FIXMATH_NO_ROUNDING
    quotient++;
    #endif
    
    fix16_t result = quotient >> 1;
    
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

#if 0
fix16_t fix16_div(fix16_t a, fix16_t b)
{
    // This uses the basic binary non-restoring division algorithm.
    // It appears to be faster to do the whole division manually than
    // trying to compose a 64-bit divide out of 32-bit divisions and
    // multiplications, especially on platforms without hardware divide.
    
    if (b == 0)
        return fix16_min;
    
    int32_t remainder = a;
    uint32_t divider = (b >= 0) ? b : (-b);
    uint32_t p_quotient = 0; // positive bits
    uint32_t n_quotient = 0; // negative bits
    uint32_t bit_pos = 65536;
    
    // The non-restoring division requires that |D| >= |N|.
    uint32_t abs_rem = (remainder >= 0) ? remainder : (-remainder);
    while (divider < abs_rem)
    {
        divider <<= 1;
        bit_pos <<= 1;
    }
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (bit_pos == 0)
        return fix16_overflow;
    #endif
    bit_pos >>= 1;
    
    if (divider & 0x80000000)
    {
        // Divider would not fit in signed 32 bit integer.
        // Perform one step of division separately to avoid overflows later.
        // We know that the dividers bottom bit is zero here.
        divider >>= 1;
        if (remainder > 0)
        {
            p_quotient |= bit_pos;
            remainder -= divider;
        }
        else
        {
            n_quotient |= bit_pos;
            remainder += divider;
        }
        bit_pos >>= 1;
    }
    
    // This is the division main loop
    while (bit_pos)
    {
        if (remainder > 0)
        {
            p_quotient |= bit_pos;
            
            // This shift may overflow, but it all works out after the
            // subtraction :)
            remainder <<= 1;
            remainder -= divider;
        }
        else if (remainder < 0)
        {
            n_quotient |= bit_pos;
            
            remainder <<= 1;
            remainder += divider;
        }
        else
        {
            break; // All done, remainder = 0
        }
        
        bit_pos >>= 1;
    }
    
    // Convert result to normal binary
    int32_t result = p_quotient - n_quotient;
    
    #ifndef FIXMATH_NO_ROUNDING
    // Figure out rounding based on the remainder
    // This is the same as calculating two extra bits, like this:
    // 00  -3  1/4
    // 0-  -2  2/4
    // 01  -1  3/4
    // --   0  0/4
    // 10   1  1/4
    // 1-   2  2/4
    // 11   3  3/4
    // ^    ^  ^--- fractional part of the result (2/4, 3/4 round away from 0)
    // |    `------ extra bits converted to decimal, neg -> result--
    // `----------- extra bits (+1,-1 coded)
    
    if (remainder != 0)
    {
        if (remainder > 0)
        {
            // First extra bit is +1
            remainder <<= 1;
            remainder -= divider;
        }
        else if (remainder < 0)
        {
            // First extra bit is -1
            result--;
            remainder <<= 1;
            remainder += divider;
        }
        
        if (result >= 0 && remainder >= 0)
            result++;
        else if (result < 0 && remainder > 0)
            result++;
    }
    #endif
    
    if (b < 0)
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




