#include "fix16_base.h"

fix16_t fix16_omul(fix16_t a, fix16_t b)
{
    int64_t product = (int64_t)a * b;
    
    // The upper 17 bits should all be the same (the sign).
    if (product >> 63 != product >> 47)
        return fix16_overflow;
    
#ifdef FIXMATH_NO_ROUNDING
    return product >> 16;
#else
    // Subtracting 0x8000 (= 0.5) and then using signed right shift
    // achieves proper rounding to result-1, except in the corner
    // case of negative numbers and lowest word = 0x8000.
    // To handle that, we also have to subtract 1 for negative numbers.
    product -= 0x8000;
    product -= (uint64_t)product >> 63;
    
    // Discard the lowest 16 bits. Note that this is not exactly the same
    // as dividing by 0x10000. For example if product = -1, result will
    // also be -1 and not 0. This is compensated by adding +1 to the result
    // and compensating this in turn in the rounding above.
    fix16_t result = product >> 16;
    result += 1;
    return result;
#endif
}

fix16_t fix16_odiv(fix16_t a, fix16_t b)
{
#ifdef FIXMATH_NO_ROUNDING
    int64_t quotient = ((int64_t)a << 16) / b;
    if (quotient >> 63 != quotient >> 31) return fix16_overflow;
    return quotient;
#else
    // To implement rounding, we first shift temp left by 17 bits instead of
    // the 16 required by the fixed point format. The lowest bit is used
    // for rounding.
    // a/b +- 0.5 = (2a/b +- 1)/2
    int64_t a_x17 = (int64_t)a << 17;
    int64_t quotient = a_x17 / b;
    
    if (quotient >> 63 != quotient >> 32)
        return fix16_overflow;
    
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




