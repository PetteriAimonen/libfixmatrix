#include "fixmatrix.h"
#include "fix16_base.h"

// Helper macros for detecting overflows in basic operations.
// To reduce the number of conditional branches, they use bitwise
// logic instead of if's. The overflow status of all computations
// is accumulated in a variable and checked at the end of a loop.

#define SIGN(x) ((uint32_t)(x) >> 31)
#define NOT(x) ((x) ^ 1)

// Detect overflows in multiplication
// Sign of product should be sign of a ^ sign of b
// Note: Works only for product != 0
#define MUL_OVERFLOW(a, b, product) (SIGN(a)^SIGN(b)^SIGN(product))

// Detect overflows in division
// Sign of result should be sign of divider ^ sign of value
// Note: works only for result != 0
#define DIV_OVERFLOW(a, b, result) (SIGN(a)^SIGN(b)^SIGN(result))

// Detect overflows in addition
// overflow can only happen if sign of a == sign of b,
// and then it causes sign of sum != sign of a            
#define ADD_OVERFLOW(a, b, sum) (NOT(SIGN(a) ^ SIGN(b))                           & (SIGN(sum) ^ SIGN(a)))

// Detect overflows in subtraction
// Overflow can occur if sign of a != sign of b
// and is detected by sign of diff != sign of a
#define SUB_OVERFLOW(a, b, diff) ((SIGN(a) ^ SIGN(b)) & (SIGN(a) ^ SIGN(diff)))


/*********************************
 * Operations between 2 matrices *
 *********************************/

// Calculates the dotproduct of two vectors of size n.
// If overflow happens, sets flag in errors
fix16_t dotproduct(const fix16_t *a, uint8_t a_stride, const fix16_t *b, uint8_t b_stride, uint8_t n, uint8_t *errors)
{
    fix16_t sum = 0;
    uint8_t overflows = 0;
    
    while (n--)
    {
        // Compute result
        fix16_t product = fix16_mul(*a, *b);
        
        if (product != 0)
        {
            fix16_t newsum = sum + product;
            overflows |= MUL_OVERFLOW(*a, *b, product);
            overflows |= ADD_OVERFLOW(sum, product, newsum);
            sum = newsum;
        }
        
        // Go to next item
        a += a_stride;
        b += b_stride;
    }
    
    if (overflows)
    {
        *errors |= FIXMATRIX_OVERFLOW;
    }
    
    return sum;
}

void mf16_mul(mf16 *dest, const mf16 *a, const mf16 *b)
{
    int row, column;
    dest->errors = a->errors | b->errors;
    
    if (a->columns != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
    if (dest == a || dest == b)
        dest->errors |= FIXMATRIX_USEERR;
    
    dest->rows = a->rows;
    dest->columns = b->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = dotproduct(
                &a->data[row][0], 1,
                &b->data[0][column], FIXMATRIX_MAX_SIZE,
                a->columns, &dest->errors);
        }
    }
}

// Multiply transpose of at with b
void mf16_mul_t(mf16 *dest, const mf16 *at, const mf16 *b)
{
    int row, column;
    dest->errors = at->errors | b->errors;
    
    if (at->rows != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
    if (dest == at || dest == b)
        dest->errors |= FIXMATRIX_USEERR;
    
    dest->rows = at->columns;
    dest->columns = b->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = dotproduct(
                &at->data[0][row], FIXMATRIX_MAX_SIZE,
                &b->data[0][column], FIXMATRIX_MAX_SIZE,
                at->rows, &dest->errors);
        }
    }
}

void mf16_add(mf16 *dest, const mf16 *a, const mf16 *b)
{
    int row, column;
    uint8_t overflows = 0;
    
    dest->errors = a->errors | b->errors;
    if (a->columns != b->columns || a->rows != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
    dest->rows = a->rows;
    dest->columns = a->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            fix16_t _a = a->data[row][column];
            fix16_t _b = b->data[row][column];
            fix16_t sum = _a + _b;
            overflows |= ADD_OVERFLOW(_a, _b, sum);
            dest->data[row][column] = sum;
        }
    }
    
    if (overflows)
    {
        dest->errors |= FIXMATRIX_OVERFLOW;
    }
}

void mf16_sub(mf16 *dest, const mf16 *a, const mf16 *b)
{
    int row, column;
    uint8_t overflows = 0;
    
    dest->errors = a->errors | b->errors;
    if (a->columns != b->columns || a->rows != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
    dest->rows = a->rows;
    dest->columns = a->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            fix16_t _a = a->data[row][column];
            fix16_t _b = b->data[row][column];
            fix16_t diff = _a - _b;
            overflows |= SUB_OVERFLOW(_a, _b, diff);
            dest->data[row][column] = diff;
        }
    }
    
    if (overflows)
    {
        dest->errors |= FIXMATRIX_OVERFLOW;
    }
}

/*********************************
 * Operations on a single matrix *
 *********************************/

void mf16_transpose(mf16 *dest, const mf16 *matrix)
{
    int row, column;
    
    // This code is a bit tricky in order to work
    // in the situation when dest = matrix.
    // Before writing a value in dest, we must copy
    // the corresponding value from matrix to a temporary
    // variable.
    
    // We actually transpose a n by n square matrix, because
    // that can be done in-place easily. Because mf16 always
    // allocates a square area even if actual matrix is smaller,
    // this is not a problem.
    int n = matrix->rows;
    if (matrix->columns > n) n = matrix->columns;
    
    uint8_t rows = matrix->rows;
    dest->rows = matrix->columns;
    dest->columns = rows;
    dest->errors = matrix->errors;
    
    for (row = 0; row < n; row++)
    {
        for (column = 0; column < row; column++)
        {
            fix16_t temp = matrix->data[row][column];
            dest->data[row][column] = matrix->data[column][row];
            dest->data[column][row] = temp;
        }
        
        dest->data[row][row] = matrix->data[row][row];
    }
}

/***************************************
 * Operations of a matrix and a scalar *
 ***************************************/

void mf16_mul_s(mf16 *dest, const mf16 *matrix, fix16_t scalar)
{
    int row, column;
    uint8_t overflows = 0;
    
    dest->rows = matrix->rows;
    dest->columns = matrix->columns;
    dest->errors = matrix->errors;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            fix16_t value = matrix->data[row][column];
            fix16_t product = fix16_mul(value, scalar);
            
            if (product != 0)
            {
                overflows |= MUL_OVERFLOW(value, scalar, product);
            }
            
            dest->data[row][column] = product;
        }
    }
    
    if (overflows)
    {
        dest->errors |= FIXMATRIX_OVERFLOW;
    }
}


/***************************************************
 * Solving linear equations using QR decomposition *
 ***************************************************/

// Takes two columns vectors, v and u, of size n.
// Performs v = v - dot(u, v) * u,
// where dot(u,v) has already been computed
// u is assumed to be an unit vector.
static void subtract_projection(fix16_t *v, const fix16_t *u, fix16_t dot, int n, uint8_t *errors)
{
    uint8_t overflows = 0;
    
    while (n--)
    {
        // For unit vector u, u[i] <= 1
        // Therefore this multiplication cannot overflow
        fix16_t product = fix16_mul(dot, *u);
        
        // TODO: Can this overflow ever happen?
        fix16_t diff = *v - product;
        overflows |= SUB_OVERFLOW(*v, product, diff);
        
        *v = diff;
        
        v += FIXMATRIX_MAX_SIZE;
        u += FIXMATRIX_MAX_SIZE;
    }
    
    if (overflows)
    {
        *errors |= FIXMATRIX_OVERFLOW;
    }
}

void mf16_qr_decomposition(mf16 *q, mf16 *r, const mf16 *matrix, int reorthogonalize)
{
    int i, j, reorth;
    fix16_t dot, norm;
    
    uint8_t stride = FIXMATRIX_MAX_SIZE;
    uint8_t n = matrix->rows;
    
    // This uses the modified Gram-Schmidt algorithm.
    // subtract_projection takes advantage of the fact that
    // previous columns have already been normalized.
    
    // We start with q = matrix
    if (q != matrix)
    {
        *q = *matrix;
    }
    
    // R is initialized to have square size of cols(A) and zeroed.
    r->columns = matrix->columns;
    r->rows = matrix->columns;
    r->errors = 0;
    for (j = 0; j < r->rows; j++)
    {
        for (i = 0; i < r->columns; i++)
        {
            r->data[j][i] = 0;
        }
    }
    
    // Now do the actual Gram-Schmidt for the rows.
    for (j = 0; j < q->columns; j++)
    {
        for (reorth = 0; reorth <= reorthogonalize; reorth++)
        {
            for (i = 0; i < j; i++)
            {
                fix16_t *v = &q->data[0][j];
                fix16_t *u = &q->data[0][i];
                
                dot = dotproduct(v, stride, u, stride, n, &q->errors);
                subtract_projection(v, u, dot, n, &q->errors);
                
                r->data[i][j] += dot;
            }
        }
        
        // Normalize the row in q
        dot = dotproduct(&q->data[0][j], stride, &q->data[0][j], stride,
                         n, &q->errors);
        norm = fix16_sqrt(dot);
        r->data[j][j] = norm;
        
        if (norm < 5 && norm > -5)
        {
            // Nearly zero norm, which means that the row
            // was linearly dependent.
            q->errors |= FIXMATRIX_SINGULAR;
            continue;
        }
        
        for (i = 0; i < n; i++)
        {
            // norm >= v[i] for all i, therefore this division
            // doesn't overflow unless norm approaches 0.
            q->data[i][j] = fix16_div(q->data[i][j], norm);
        }
    }
    
    r->errors = q->errors;
}

void mf16_solve(mf16 *dest, const mf16 *q, const mf16 *r, const mf16 *matrix)
{
    int row, column, variable;
    uint8_t overflows = 0;
    
    // Ax=b <=> QRx=b <=> Q'QRx=Q'b <=> Rx=Q'b
    // Q'b is calculated directly and x is then solved row-by-row.
    mf16_mul_t(dest, q, matrix);
    
    if (r->columns != r->rows || r->columns != q->columns)
    {
        dest->errors |= FIXMATRIX_USEERR;
    }
    
    for (column = 0; column < dest->columns; column++)
    {
        for (row = dest->rows - 1; row >= 0; row--)
        {
            fix16_t value = dest->data[row][column];
            
            // Subtract any already solved variables
            for (variable = row + 1; variable < r->columns; variable++)
            {
                fix16_t multiplier = r->data[row][variable];
                fix16_t known_value = dest->data[variable][column];
                fix16_t product = fix16_mul(multiplier, known_value);
                
                if (product != 0)
                {
                    fix16_t newvalue = value - product;
                
                    overflows |= MUL_OVERFLOW(multiplier, known_value, product);
                    overflows |= SUB_OVERFLOW(value, product, newvalue);
                    
                    value = newvalue;
                }
            }
            
            // Now value = R_ij x_i <=> x_i = value / R_ij
            fix16_t divider = r->data[row][row];
            if (divider == 0)
            {
                dest->errors |= FIXMATRIX_SINGULAR;
                continue;
            }
            
            fix16_t result = fix16_div(value, divider);
            dest->data[row][column] = result;
            
            if (result != 0)
            {
                overflows |= DIV_OVERFLOW(value, divider, result);
            }
        }
    }
    
    if (overflows)
    {
        dest->errors |= FIXMATRIX_OVERFLOW;
    }
}


