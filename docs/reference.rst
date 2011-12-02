================================
Libfixmatrix: Function reference
================================

.. contents ::

fixmatrix.h
===========

FIXMATRIX_MAX_SIZE
------------------
Maximum size of matrices. Configurable option. ::

    #ifndef FIXMATRIX_MAX_SIZE
    #define FIXMATRIX_MAX_SIZE 8
    #endif

This option defines the allocation size of matrices. The matrices and vectors
can be smaller than this and also non-square, but a memory is always allocated for a matrix of this size.

For the default value of 8, each matrix takes 260 bytes of RAM.

You can modify the option by editing the header or by setting it on compiler command line.

mf16
----
Data structure for matrices. ::

    typedef struct {
        uint8_t rows;
        uint8_t columns;
        uint8_t errors;
        fix16_t data[FIXMATRIX_MAX_SIZE][FIXMATRIX_MAX_SIZE];
    } mf16;

:rows:      Number of rows in the matrix, 1 <= rows <= FIXMATRIX_MAX_SIZE.
:columns:   Number of columns in the matrix, 1 <= columns <= FIXMATRIX_MAX_SIZE.
:errors:    Error flags for the result. Should be initialized to 0 when matrix is filled with data, and checked for != 0 in output. Any flags that are set will be propagated through operations.
:data:      Numeric values in the matrix. Stored in row-major format, i.e. data[row][column].

A matrix can be created and initialized with data using standard C syntax::

    mf16 mymatrix = {2, 2, 0,
        {{1, 0}, // First row
         (0, 1)} // Second row
    };

.. sidebar:: Handling of overflows
    
    The *fix16_t* datatype has a value range from -32768.0 to +32767.0. Any
    results outside this range will be corrupted, but this is automatically
    detected.
    
    The happening of an overflow is indicated by setting *FIXMATRIX_OVEFLOW*
    in the result matrix. You typically want to discard the result if this
    happens.
    
    Most functions in libfixmatrix overflow only if the result exceeds the
    maximum range. However, in `mf16_qr_decomposition`_ and `mf16_solve`_
    overflows can happen even if the result would fit the datatype.
    This is explained in more detail in the documentation of these functions.
    
Error flags
-----------
The following error flags have been defined for mf16.errors::

    #define FIXMATRIX_OVERFLOW 0x01
    #define FIXMATRIX_DIMERR   0x02
    #define FIXMATRIX_USEERR   0x04
    #define FIXMATRIX_SINGULAR 0x08

Each flag is set when a particular condition occurs. The flag is set in the
result matrix, while the operands are left unmodified (unless, of course,
the operand is the same as the destination).

:OVERFLOW:  An intermediate result has exceeded the range of *fix16_t* data type, i.e. -32768.0 to +32767.0.
:DIMERR:    Matrices passed in have incompatible dimensions. E.g. addition requires that the matrices are of same size, and multiplication requires that a.cols == b.rows.
:USEERR:    The function is used in an unsupported way. E.g. matrices alias in multiplication or the QR decomposition passed to `mf16_solve`_ is not actually a QR composition.
:SINGULAR:  In `mf16_solve`_ and `mf16_qr_decomposition`_, the matrix was rank-deficient and therefore not invertible.

.. sidebar:: Pointer aliasing

    Most of the functions in libfixmatrix take a pointer for both the
    destination, where result will be stored, and the operands of the
    computation. Pointer aliasing means that the destination actually
    points to one of the operands, i.e. the equivalent of a = a + b.
    
    In most cases this is allowed and is a good practice, because it
    reduces the amount of matrices you need to allocate. However,
    multiplication and `mf16_solve`_ don't support this.
    
    The read-only (const \*) operands can always alias with other read-only
    operands. Aliasing only matters for the destination operands.
    
    The function descriptions point out whether aliasing is allowed or not.
    Additionally, each function where it is not allowed checks for improper
    aliasing and sets *FIXMATRIX_USEERR* if you call it in a wrong way.

dotproduct
----------
Calculates the dot product of two sequences of *fix16_t* numbers::

    fix16_t dotproduct(const fix16_t *a, uint8_t a_stride,
                       const fix16_t *b, uint8_t b_stride,
                       uint8_t n, uint8_t *errors);

:a:         Pointer to the first number of the first sequence.
:a_stride:  Increment to the next number of the sequence, specified in terms of *sizeof(fix16_t)*. I.e. \*(a + a_stride) is the second entry in first sequence.
:b:         Second sequence.
:b_stride:  Stride of the second sequence.
:n:         Number of entries in each sequence.
:errors:    Pointer to variable where *FIXMATRIX_OVERFLOW* will be set if the result overflows.
:returns:   The dot product of a and b, that is, each entry of a multiplied by the corresponding entry of b and summed together.

mf16_mul
--------
Matrix multiplication, dest = a * b::
    
    void mf16_mul(mf16 *dest, const mf16 *a, const mf16 *b);

:dest:      Destination for storing the result. Cannot alias with *a* or *b*.
:a:         Left operand of the multiplication.
:b:         Right operand of the multiplication.

Matrix multiplication requires that the number of rows in *b* equals the number of columns in *a*. If this is not the case, FIXMATRIX_DIMERR is set.

Result will have *a->rows* rows and *b->columns* columns.

mf16_mul_t
----------
Matrix multiplication where the first argument is transposed, dest = a' * b::

    void mf16_mul_t(mf16 *dest, const mf16 *at, const mf16 *b);

:dest:      Destination for storing the result. Cannot alias with *at* or *b*.
:at:        Left operand of the multiplication. Will be used in a transposed order.
:b:         Right operand of the multiplication.

The number of rows in *b* must equal the number of rows in *at*.
Result will have *at->columns* rows and *b->columns* columns.

mf16_add
--------
Matrix addition, dest = a + b::

    void mf16_add(mf16 *dest, const mf16 *a, const mf16 *b);
    
:dest:      Destination for storing the result. Can be same as *a* or *b* or both.
:a:         First matrix in addition.
:b:         Second matrix in addition.

The matrices are added entry-by-entry. The matrices *a* and *b* must have the same dimensions.

mf16_sub
--------
Matrix subtraction, dest = a - b:

    void mf16_sub(mf16 *dest, const mf16 *a, const mf16 *b);

:dest:      Destination for storing the result. Can be same as *a* or *b* or both.
:a:         Matrix to subtract from.
:b:         Matrix to subtract.

Each entry of *b* is subtracted from the corresponding entry in *a*. Matrices
must have the same dimensions.

mf16_transpose
--------------
Transposition of a matrix, dest = matrix'::

    void mf16_transpose(mf16 *dest, const mf16 *matrix);

:dest:      Destination for storing the result. Can be same as *matrix*.
:matrix:    Matrix to transpose. Can have any dimensions.

mf16_mul_s
----------
Multiplication of matrix by scalar, dest = s * matrix::

    void mf16_mul_s(mf16 *dest, const mf16 *matrix, fix16_t scalar);

:dest:      Destination for storing the result. Can be same as *matrix*.
:matrix:    Matrix to multiply.
:scalar:    Scalar value to multiply by.

Each entry of *matrix* is multiplied by the scalar value.
    
mf16_qr_decomposition
---------------------
QR-decomposition of a matrix, q * r = matrix::

    void mf16_qr_decomposition(mf16 *q, mf16 *r, const mf16 *matrix, int reorthogonalize);

:q:         Destination for the orthonormal part of the result.
:r:         Destination for the upper-triangular part of the result.
:matrix:    Matrix to decompose.
:reorthogonalize: Iteration count, larger values improve precision. Value of 0 is fastest and gives usually error of less than 0.1%. If rounding is not disabled (by defining *FIXMATH_NO_ROUNDING*), values larger than 1 don't improve precision. If rounding is disabled, values up to 3 may be useful.

QR-decomposition is the first phase of solving an equation system using libfixmatrix.
It can be used both for exact solutions using square matrices and for least squares solutions with rectangular matrices.

One of the destination matrices *q* and *r* may alias with *matrix*. The execution
time is the shortest when *q* = *matrix*, because that avoids one matrix-sized memory copy inside the function.

*Matrix* should have a rank equal to the number of columns, i.e. it should have a full column rank, i.e. all of its columns should be linearly independent. A matrix that does not have full column rank does not have an unique solution. This function will report that by setting *FIXMATRIX_SINGULAR* in both of the result matrices.

When *matrix* contains more rows than columns, an economy factorization is returned.
This means that q is not a square matrix, but otherwise the usual properties of QR-decomposition hold.

The values in the *q* matrix are, by definition, less than or equal to 1.
Therefore they can have at most 16 bits of precision. This naturally limits
the precision obtained from any further calculations, which may be important
if the matrix in question contains large values. It would be possible to use
a different fixed-point scaling for the Q matrix, but it would increase code
size and is not currently implemented.

Note that this function will cause overflows before the values in matrices
involved are even close to the *fix16_t* limits. This is because it internally
computes the norm of each column, and the intermediate product of this is the
square of the norm. Therefore overflows may happen if the norm of any column
exceeds 256.0. The overflow may not happen every time, though, because the
linearly dependent parts of the column are subtracted first.

