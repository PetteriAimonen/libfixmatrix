/* Basic quaternion implementation based on libfixmath fix16_t datatype. */

#ifndef _FIXQUAT_H_
#define _FIXQUAT_H_

#include <fix16.h>
#include "fixmatrix.h"

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

// Multiply quaternion by scalar
void qf16_mul_s(qf16 *dest, const qf16 *q, fix16_t s);

// Divide quaternion by scalar
void qf16_div_s(qf16 *dest, const qf16 *q, fix16_t s);

// Quaternion norm
fix16_t qf16_norm(const qf16 *q);

// Normalize quaternion
void qf16_normalize(qf16 *dest, const qf16 *q);

// Unit quaternion to rotation matrix
void qf16_to_matrix(mf16 *dest, const qf16 *q);

#endif