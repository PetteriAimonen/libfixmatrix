#ifndef _FIX16_BASE_H_
#define _FIX16_BASE_H_

#include <stdint.h>

typedef int32_t fix16_t;
static const fix16_t fix16_one = 0x10000;
static const fix16_t fix16_max = 0x7FFFFFFF;
static const fix16_t fix16_min = 0x80000000;

/* Multiplication and division */
fix16_t fix16_mul(fix16_t a, fix16_t b);
fix16_t fix16_div(fix16_t a, fix16_t b);

/* Conversion functions between fix16_t and float/integer. */
static inline fix16_t fix16_from_int(int a) { return a * fix16_one; }
static inline float fix16_to_float(fix16_t a) { return (float)a / fix16_one; }
static inline double fix16_to_double(fix16_t a) { return (double)a / fix16_one; }

static inline int fix16_to_int(fix16_t a)
{
    return (int)fix16_div(a, fix16_one);
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
