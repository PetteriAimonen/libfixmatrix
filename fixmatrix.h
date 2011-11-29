/* A very basic and small matrix algebra library atop libfixmath
 * fixed point numbers. Suitable for small matrices, usually less
 * than 10x10.
 * 
 * This library does not do saturating arithmetic, but it does
 * feature an overflow flag for detecting erroneous results.
 * 
 * Goals of the library are small footprint and fast execution.
 * Only very basic operations are supported.
 * 
 * Matrices can have any size from 1x1 up to FIXMATRIX_MAX_SIZE
 * (configurable), also non-square, but memory is always allocated
 * for the maximum size.
 * 
 * Error handling is done using flags in the matrix structure.
 * This makes it easy to detect if any errors occurred in any of
 * the computations, without checking a return status from each
 * function.
 * 
 * The functions still perform the calculations even if the result
 * is known to be erroneous.
 */

#ifndef __fixmatrix_h_
#define __fixmatrix_h_

#include <stdint.h>
#include <fix16.h>

// Maximum size of matrices.
#ifndef FIXMATRIX_MAX_SIZE
#define FIXMATRIX_MAX_SIZE 8
#endif

typedef struct {
    uint8_t rows;
    uint8_t columns;
    
    /* Error flags are used to detect exceptions in the computations.
     * The flags are automatically propagated to the result if either
     * of the operands is invalid.
     * Currently the following flags are defined:
     * - FIXMATRIX_OVERFLOW: A value has exceeded 32767 and wrapped around
     * - FIXMATRIX_DIMERR: Operands have incompatible dimensions
     * - FIXMATRIX_USEERR: Function was called in unsupported way
     */
    uint8_t errors;
    
    /* Data is stored in memory in row-major format, e.g.
     * entry at (row, column) is data[row][column]
     */
    fix16_t data[FIXMATRIX_MAX_SIZE][FIXMATRIX_MAX_SIZE];
} mf16;

#define FIXMATRIX_OVERFLOW 0x01
#define FIXMATRIX_DIMERR   0x02
#define FIXMATRIX_USEERR   0x04
#define FIXMATRIX_SINGULAR 0x08

// Calculates the dotproduct of two vectors of size n.
// If overflow happens, sets flag in errors
fix16_t dotproduct(const fix16_t *a, uint8_t a_stride,
                   const fix16_t *b, uint8_t b_stride,
                   uint8_t n, uint8_t *errors);

// Operations between two matrices
// Note: in multiplication, dest can not point to the
// same matrix as a or b. If it does, FIXMATRIX_USEERR
// is set.
void mf16_mul(mf16 *dest, const mf16 *a, const mf16 *b);
void mf16_mul_t(mf16 *dest, const mf16 *at, const mf16 *b);

// In addition and subtraction, a = dest and b = dest are allowed.
void mf16_add(mf16 *dest, const mf16 *a, const mf16 *b);
void mf16_sub(mf16 *dest, const mf16 *a, const mf16 *b);

// Operations on a single matrix
// matrix and dest can alias.
void mf16_transpose(mf16 *dest, const mf16 *matrix);

// Operations of a matrix and a scalar
// matrix and dest can alias.
void mf16_mul_s(mf16 *dest, const mf16 *matrix, fix16_t scalar);

// QR-decomposition of a matrix
//
// This function does not support rank-deficient matrices.
// If rank(A) < cols(A), FIXMATRIX_SINGULAR is set.
//
// Overdetermined systems of rows(A) > cols(A) are supported.
// For them, an 'economy' factorization is returned, with q
// being non-square. mf16_solve will then return least squares
// solution.
//
// Specifying reorthogonalize > 0 increases iterations and
// improves result accuracy. Reorthogonalize = 0 is the fastest
// and typically gives rounding error of 0.1-0.5%. Values >3
// rarely improve precision.
//
// In mf16_qr_decomposition, q and matrix, or alternatively, r and matrix
// may alias i.e. point to the same memory location.
void mf16_qr_decomposition(mf16 *q, mf16 *r, const mf16 *matrix, int reorthogonalize);

// Solving a system of linear equations Ax = b, or equivalently,
// right division A\b, using QR-decomposition.
// matrix is the b and x is stored to dest.
// The arguments cannot alias.
// matrix may have multiple columns, which are then solved
// independently.
// If you really really want and think that it is a
// good idea to invert matrices, you can do it by
// passing identity matrix as 'matrix'.
void mf16_solve(mf16 *dest, const mf16 *q, const mf16 *r, const mf16 *matrix);

#endif
