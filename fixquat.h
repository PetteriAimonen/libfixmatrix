/* Basic quaternion implementation based on libfixmath fix16_t datatype. */

#ifndef _FIXQUAT_H_
#define _FIXQUAT_H_

#include <fix16.h>
#include "fixmatrix.h"
#include "fixvector3d.h"

typedef struct {
    fix16_t a; // Real part
    fix16_t b; // i
    fix16_t c; // j
    fix16_t d; // k
} qf16;

// Conjugate of quaternion
void qf16_conj(qf16 *dest, const qf16 *q);

// Multiply two quaternions, dest = q * r.
void qf16_mul(qf16 *dest, const qf16 *q, const qf16 *r);

// Add two quaternions, dest = q + r
void qf16_add(qf16 *dest, const qf16 *q, const qf16 *r);

// Multiply quaternion by scalar
void qf16_mul_s(qf16 *dest, const qf16 *q, fix16_t s);

// Divide quaternion by scalar
void qf16_div_s(qf16 *dest, const qf16 *q, fix16_t s);

// Dot product of two quaternions
fix16_t qf16_dot(const qf16 *q, const qf16 *r);

// Quaternion norm
fix16_t qf16_norm(const qf16 *q);

// Normalize quaternion
void qf16_normalize(qf16 *dest, const qf16 *q);

// Quaternion power (exponentation)
void qf16_pow(qf16 *dest, const qf16 *q, fix16_t power);

// Weighted average of two quaternions
// Think of it as q = w * q1 + (1 - w) * q2, but the internal algorithm considers attitudes.
void qf16_avg(qf16 *dest, const qf16 *q1, const qf16 *q2, fix16_t weight);

// Unit quaternion from axis and angle.
// Axis should have unit length and angle in radians.
void qf16_from_axis_angle(qf16 *dest, const v3d *axis, fix16_t angle);

// Unit quaternion to rotation matrix
void qf16_to_matrix(mf16 *dest, const qf16 *q);

// Rotate vector using quaternion
void qf16_rotate(v3d *dest, const qf16 *q, const v3d *v);

static inline void qf16_from_v3d(qf16 *q, const v3d *v, fix16_t a)
{
    q->a = a;
    q->b = v->x;
    q->c = v->y;
    q->d = v->z;
}

static inline void qf16_to_v3d(v3d *v, const qf16 *q)
{
    v->x = q->b;
    v->y = q->c;
    v->z = q->d;
}

#endif