#include "fixvector3d.h"
#include "fixarray.h"

void v3d_add(v3d *dest, const v3d *a, const v3d *b)
{
        dest->x = fix16_add(a->x, b->x);
        dest->y = fix16_add(a->y, b->y);
        dest->z = fix16_add(a->z, b->z);
}

void v3d_sub(v3d *dest, const v3d *a, const v3d *b)
{
        dest->x = fix16_sub(a->x, b->x);
        dest->y = fix16_sub(a->y, b->y);
        dest->z = fix16_sub(a->z, b->z);
}

void v3d_mul_s(v3d *dest, const v3d *a, fix16_t b)
{
        dest->x = fix16_mul(a->x, b);
        dest->y = fix16_mul(a->y, b);
        dest->z = fix16_mul(a->z, b);
}

void v3d_div_s(v3d *dest, const v3d *a, fix16_t b)
{
        dest->x = fix16_div(a->x, b);
        dest->y = fix16_div(a->y, b);
        dest->z = fix16_div(a->z, b);
}

// Norm
fix16_t v3d_norm(const v3d *a)
{
    return fa16_norm(&a->x, &a->y - &a->x, 3);
}

void v3d_normalize(v3d *dest, const v3d *a)
{
    v3d_div_s(dest, a, v3d_norm(a));
}

// Dot product
fix16_t v3d_dot(const v3d *a, const v3d *b)
{
    return fa16_dot(&a->x, &a->y - &a->x, &b->x, &b->y - &b->x, 3);
}

// Cross product
void v3d_cross(v3d *dest, const v3d *a, const v3d *b)
{
    v3d tmp;
    fa16_unalias(dest, (void**)&a, (void**)&b, &tmp, sizeof(tmp));
    
    dest->x = fix16_sub(fix16_mul(a->y, b->z), fix16_mul(a->z, b->y));
    dest->y = fix16_sub(fix16_mul(a->z, b->x), fix16_mul(a->x, b->z));
    dest->z = fix16_sub(fix16_mul(a->x, b->y), fix16_mul(a->y, b->x));
}