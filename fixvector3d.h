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
v3d v3d_add(v3d a, v3d b);
v3d v3d_sub(v3d a, v3d b);
v3d v3d_mul(v3d a, fix16_t b);
v3d v3d_div(v3d a, fix16_t b);

// Norm
fix16_t v3d_norm(v3d a);

// Dot product
fix16_t v3d_dot(v3d a, v3d b);

// Cross product
v3d v3d_cross(v3d a, v3d b);

#endif