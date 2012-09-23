#include "fixvector2d.h"
#include "fixarray.h"

// Basic arithmetic
void v2d_add(v2d *dest, const v2d *a, const v2d *b)
{
    dest->x = fix16_add(a->x, b->x);
    dest->y = fix16_add(a->y, b->y);
}

void v2d_sub(v2d *dest, const v2d *a, const v2d *b)
{
    dest->x = fix16_sub(a->x, b->x);
    dest->y = fix16_sub(a->y, b->y);
}

void v2d_mul_s(v2d *dest, const v2d *a, fix16_t b)
{
    dest->x = fix16_mul(a->x, b);
    dest->y = fix16_mul(a->y, b);
}

void v2d_div_s(v2d *dest, const v2d *a, fix16_t b)
{
    dest->x = fix16_div(a->x, b);
    dest->y = fix16_div(a->y, b);
}

// Norm
fix16_t v2d_norm(const v2d *a)
{
    return fa16_norm(&a->x, &a->y - &a->x, 2);
}

void v2d_normalize(v2d *dest, const v2d *a)
{
    v2d_div_s(dest, a, v2d_norm(a));
}

// Dot product
fix16_t v2d_dot(const v2d *a, const v2d *b)
{
    return fix16_add(fix16_mul(a->x, b->x), fix16_mul(a->y, b->y));
}

// Rotation (positive direction = counter-clockwise, angle in radians)
void v2d_rotate(v2d *dest, const v2d *a, fix16_t angle)
{
    fix16_t c = fix16_cos(angle);
    fix16_t s = fix16_sin(angle);
    
    dest->x = fix16_add(fix16_mul(c, a->x), fix16_mul(-s, a->y));
    dest->y = fix16_add(fix16_mul(s, a->x), fix16_mul(c, a->y));
}
