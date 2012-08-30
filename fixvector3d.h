/* 3D vector operations */

#ifndef _fixvector3d_h_
#define _fixvector3d_h_

#include <fix16.h>

typedef struct {
	fix16_t x;
	fix16_t y;
	fix16_t z;
} v3d;

// Basic arithmetic
void v3d_add(v3d *dest, const v3d *a, const v3d *b);
void v3d_sub(v3d *dest, const v3d *a, const v3d *b);
void v3d_mul_s(v3d *dest, const v3d *a, fix16_t b);
void v3d_div_s(v3d *dest, const v3d *a, fix16_t b);

// Norm
fix16_t v3d_norm(const v3d *a);
void v3d_normalize(v3d *dest, const v3d *a);

// Dot product
fix16_t v3d_dot(const v3d *a, const v3d *b);

// Cross product
void v3d_cross(v3d *dest, const v3d *a, const v3d *b);

#endif