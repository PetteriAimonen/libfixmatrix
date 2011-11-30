#include "fix16_base.h"

fix16_t fix16_mul(fix16_t a, fix16_t b)
{
    int64_t product = (int64_t)a * b;
    
#ifdef FIXMATH_NO_ROUNDING
    return product / 0x10000;
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

fix16_t fix16_div(fix16_t a, fix16_t b)
{
    // To implement rounding, we first shift temp left by 17 bits instead of
    // the 16 required by the fixed point format. The lowest bit is used
    // for rounding.
    // a/b +- 0.5 = (2a/b +- 1)/2
    int64_t a_x17 = (int64_t)a << 17;
    int64_t quotient = a_x17 / b;
    
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
}
