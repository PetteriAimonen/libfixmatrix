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
#include <stdbool.h>
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
#define FIXMATRIX_NEGATIVE 0x10

// Initialization functions. These expect rows and column counts to be set be the caller,
// everything else is initialized by the functions.

// Fill all the entries with the same value, and clear error status.
void mf16_fill(mf16 *dest, fix16_t value);

// Fill the diagonal entries with the given value and everything else with zeroes, and clear error status.
void mf16_fill_diagonal(mf16 *dest, fix16_t value);

// Operations between two matrices
void mf16_mul(mf16 *dest, const mf16 *a, const mf16 *b);

// Multiply transpose of at with b
void mf16_mul_at(mf16 *dest, const mf16 *at, const mf16 *b);

// Multiply a with transpose of bt
void mf16_mul_bt(mf16 *dest, const mf16 *a, const mf16 *bt);

// In addition and subtraction, a = dest and b = dest are allowed.
void mf16_add(mf16 *dest, const mf16 *a, const mf16 *b);
void mf16_sub(mf16 *dest, const mf16 *a, const mf16 *b);

// Operations on a single matrix
// matrix and dest can alias.
void mf16_transpose(mf16 *dest, const mf16 *matrix);

// Operations of a matrix and a scalar
// matrix and dest can alias.
void mf16_mul_s(mf16 *dest, const mf16 *matrix, fix16_t scalar);
void mf16_div_s(mf16 *dest, const mf16 *matrix, fix16_t scalar);

// QR-decomposition of a matrix
//
// Finds Q and R so that QR = A and Q is orthogonal and R is upper triangular.
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
// and typically gives rounding error of 0.1-0.5%. Values >1
// rarely improve precision.
//
// In mf16_qr_decomposition, q and matrix, or alternatively, r and matrix
// may alias i.e. point to the same memory location.
void mf16_qr_decomposition(mf16 *q, mf16 *r, const mf16 *matrix, int reorthogonalize);

// Solving a system of linear equations Ax = b, or equivalently,
// left division A\b, using QR-decomposition.
// matrix is the b and x is stored to dest.
// Dest can alias with matrix or q, but not with r.
// matrix may have multiple columns, which are then solved
// independently.
// If you really really want and think that it is a
// good idea to invert matrices, you can do it by
// passing identity matrix as 'matrix'.
void mf16_solve(mf16 *dest, const mf16 *q, const mf16 *r, const mf16 *matrix);

// Cholesky decomposition of a symmetric positive-definite matrix (matrix square root)
//
// Finds L so that L L' = A and L is lower triangular.
//
// Any negative square roots in computation are floored to 
// zero. If they are smaller than -0.001, FIXMATRIX_NEGATIVE
// error flag is set. Small negative values are often caused
// by rounding errors, while large negative values indicate
// non-positive definite matrix.
//
// Matrix is not checked for symmetricity. Only values in
// the lower left triangle are used.
//
// Dest and matrix can alias.
void mf16_cholesky(mf16 *dest, const mf16 *matrix);

// Matrix inversion of a matrix through its lower triangular decomposition.
//
// Finds inv(A) through L, so that A inv(A) = I and I is the identitiy matrix 
// and L a lower triangular matrix such that L L' = A.
//
// The required lower triangular matrix can be obtained through mf16_cholesky().
//
// Matrix is not checked for symmetricity. Only values in
// the lower left triangle are used.
//
// Dest and matrix can alias.
void mf16_invert_lt(mf16 *dest, const mf16 *matrix);

#endif
