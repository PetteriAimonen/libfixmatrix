/* Generic subroutines for processing fix16_t arrays (vectors). */

#ifndef _FIXARRAY_H_
#define _FIXARRAY_H_

#include <fix16.h>

// Calculates the dotproduct of two vectors of size n.
// If overflow happens, returns fix16_overflow.
fix16_t fa16_dot(const fix16_t *a, uint_fast8_t a_stride,
                 const fix16_t *b, uint_fast8_t b_stride,
                 uint_fast8_t n);

// Calculates the norm of a vector of size n.
fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n);

// Unalias function arguments using a temporary storage if necessary
// (not really related to arrays, but common to fixquat/fixvector/fixmatrix)
void fa16_unalias(void *dest, void **a, void **b, void *tmp, unsigned size);

#endif
