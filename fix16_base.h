#ifndef _FIX16_BASE_H_
#define _FIX16_BASE_H_

#include <stdint.h>

typedef int32_t fix16_t;
static const fix16_t fix16_one = 0x10000;
static const fix16_t fix16_max = 0x7FFFFFFF;
static const fix16_t fix16_min = 0x80000000;
static const fix16_t fix16_overflow = 0x80000000;

/* Multiplication and division with overflow detection.
 * These will return fix16_overflow if the result overflows.
 */
fix16_t fix16_omul(fix16_t a, fix16_t b);
fix16_t fix16_odiv(fix16_t a, fix16_t b);

/* Square root */
fix16_t fix16_sqrt(fix16_t a);

/* Subtraction and addition with overflow detection. */
static inline fix16_t fix16_oadd(fix16_t a, fix16_t b)
{
    // Use unsigned integers because overflow with signed integers is
    // an undefined operation.
    uint32_t _a = a, _b = b;
    uint32_t sum = _a + _b;
    
    // Overflow can only happen if sign of a == sign of b, and then
    // it causes sign of sum != sign of a.
    if (!((_a ^ _b) & 0x80000000) && ((_a ^ sum) & 0x80000000))
        return fix16_overflow;
    
    return sum;
}

static inline fix16_t fix16_osub(fix16_t a, fix16_t b)
{
    // Cannot invert fix16_min because of 2's complement limit asymmetry.
    if (b == fix16_min)
        return fix16_overflow;
    
    return fix16_oadd(a, -b);
}

/* Conversion functions between fix16_t and float/integer. */
static inline fix16_t fix16_from_int(int a) { return a * fix16_one; }
static inline float fix16_to_float(fix16_t a) { return (float)a / fix16_one; }
static inline double fix16_to_double(fix16_t a) { return (double)a / fix16_one; }

static inline int fix16_to_int(fix16_t a)
{
    if (a >= 0)
        return (a + fix16_one / 2) / fix16_one;
    else
        return (a - fix16_one / 2) / fix16_one;
}

static inline fix16_t fix16_from_float(float a)
{
    float temp = a * fix16_one;
#ifndef FIXMATH_NO_ROUNDING
    temp += (temp >= 0) ? 0.5f : -0.5f;
#endif
    return (fix16_t)temp;
}

static inline fix16_t fix16_from_double(double a)
{
    double temp = a * fix16_one;
#ifndef FIXMATH_NO_ROUNDING
    temp += (temp >= 0) ? 0.5f : -0.5f;
#endif
    return (fix16_t)temp;
}

#endif
