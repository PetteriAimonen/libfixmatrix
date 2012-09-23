/* 2D vector operations */

#ifndef _fixvector2d_h_
#define _fixvector2d_h_

#include <fix16.h>

typedef struct {
    fix16_t x;
    fix16_t y;
} v2d;

// Basic arithmetic
void v2d_add(v2d *dest, const v2d *a, const v2d *b);
void v2d_sub(v2d *dest, const v2d *a, const v2d *b);
void v2d_mul_s(v2d *dest, const v2d *a, fix16_t b);
void v2d_div_s(v2d *dest, const v2d *a, fix16_t b);

// Norm
fix16_t v2d_norm(const v2d *a);
void v2d_normalize(v2d *dest, const v2d *a);

// Dot product
fix16_t v2d_dot(const v2d *a, const v2d *b);

// Rotation (positive direction = counter-clockwise, angle in radians)
void v2d_rotate(v2d *dest, const v2d *a, fix16_t angle);

#endif
