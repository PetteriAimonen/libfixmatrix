#include "fixmatrix.h"

#define SIGNBIT 0x80000000

/*********************************
 * Operations between 2 matrices *
 *********************************/

// Calculates the dotproduct of row and column vectors of size n.
// If overflow happens, sets flag in errors
static fix16_t dotproduct_rowcol(const fix16_t *row, const fix16_t *column, int n, uint8_t *errors)
{
    fix16_t sum = 0;
    
    // To reduce conditional branches, the overflow status is collected
    // into this variable using bitwise operations.
    fix16_t overflows = 0;
    
    while (n--)
    {
        // Compute result
        fix16_t a = *row;
        fix16_t b = *column;
        fix16_t product = fix16_mul(a, b);
        
        if (product != 0)
        {
            fix16_t newsum = sum + product;
            
            // Detect overflows in multiplication
            // highest bit of a^b should equal the sign bit of product
            overflows |= a^b^product;
            
            // Detect overflows in addition
            // overflow can only happen if sign of product == sign of sum,
            // and then it causes sign of newsum != sign of sum
            overflows |= ~(sum ^ product) & (sum ^ newsum);
            
            sum = newsum;
        }
        
        // Go to next item
        row += 1;
        column += FIXMATRIX_MAX_SIZE;
    }
    
    if (overflows & SIGNBIT)
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
            dest->data[row][column] = dotproduct_rowcol(&a->data[row][0], &b->data[0][column], a->columns, &dest->errors);
        }
    }
}

void mf16_add(mf16 *a, const mf16 *b)
{
    int row, column;
    fix16_t overflows = 0;
    a->errors |= b->errors;
    
    if (a->columns != b->columns || a->rows != b->rows)
        a->errors |= FIXMATRIX_DIMERR;
    
    for (row = 0; row < a->rows; row++)
    {
        for (column = 0; column < a->columns; column++)
        {
            fix16_t _a = a->data[row][column];
            fix16_t _b = b->data[row][column];
            fix16_t sum = _a + _b;
            
            overflows |= ~(_a ^ _b) & (_a ^ sum);
            
            a->data[row][column] = sum;
        }
    }
    
    if (overflows & SIGNBIT)
    {
        a->errors |= FIXMATRIX_OVERFLOW;
    }
}

void mf16_sub(mf16 *a, const mf16 *b)
{
    int row, column;
    fix16_t overflows = 0;
    a->errors |= b->errors;
    
    if (a->columns != b->columns || a->rows != b->rows)
        a->errors |= FIXMATRIX_DIMERR;
    
    for (row = 0; row < a->rows; row++)
    {
        for (column = 0; column < a->columns; column++)
        {
            fix16_t _a = a->data[row][column];
            fix16_t _b = b->data[row][column];
            fix16_t diff = _a - _b;
            
            // Overflow can occur if sign of a != sign of b
            // and is detected by sign of diff != sign of a
            overflows |= (_a ^ _b) & (_a ^ diff);
            
            a->data[row][column] = diff;
        }
    }
    
    if (overflows & SIGNBIT)
    {
        a->errors |= FIXMATRIX_OVERFLOW;
    }
}

/*********************************
 * Operations on a single matrix *
 *********************************/

void mf16_transpose(mf16 *matrix)
{
    int row, column;
    
    int n = matrix->rows;
    if (matrix->columns > n) n = matrix->columns;
    
    for (row = 0; row < n; row++)
    {
        for (column = 0; column < row; column++)
        {
            fix16_t temp = matrix->data[row][column];
            matrix->data[row][column] = matrix->data[column][row];
            matrix->data[column][row] = temp;
        }
    }
    
    uint8_t rows = matrix->rows;
    matrix->rows = matrix->columns;
    matrix->columns = rows;
}

/***************************************
 * Operations of a matrix and a scalar *
 ***************************************/

void mf16_mul_s(mf16 *matrix, fix16_t scalar)
{
    int row, column;
    fix16_t overflows = 0;
    
    for (row = 0; row < matrix->rows; row++)
    {
        for (column = 0; column < matrix->columns; column++)
        {
            fix16_t value = matrix->data[row][column];
            fix16_t product = fix16_mul(value, scalar);
            
            // Detect overflows in multiplication
            // highest bit of value^scalar should equal the sign bit of product
            overflows |= value^scalar^product;
            
            matrix->data[row][column] = product;
        }
    }
    
    if (overflows & SIGNBIT)
    {
        matrix->errors |= FIXMATRIX_OVERFLOW;
    }
}


/***************************************************
 * Solving linear equations using QR decomposition *
 ***************************************************/

static fix16_t dotproduct_row(const fix16_t *v, const fix16_t *u, int n, uint8_t *errors)
{
    int i;
    fix16_t sum = 0; // Stores the dot product
    fix16_t overflows = 0;
    
    for (i = 0; i < n; i++)
    {
        fix16_t product = fix16_mul(v[i], u[i]);
        fix16_t newsum = sum + product;
        
        overflows |= v[i]^u[i]^product;
        overflows |= ~(sum ^ product) & (sum ^ newsum);
        
        sum = newsum;
    }
    
    if (overflows & SIGNBIT)
    {
        *errors |= FIXMATRIX_OVERFLOW;
    }
    
    return sum;
}

// Takes two row vectors, v and u, of size n.
// Performs v = v - dot(u, v) * u,
// where dot(u,v) has already been computed
// u is assumed to be an unit vector.
static void subtract_projection(fix16_t *v, const fix16_t *u, fix16_t dot, int n, uint8_t *errors)
{
    int i;
    fix16_t overflows = 0;
    
    for (i = 0; i < n; i++)
    {
        // For unit vector u, u[i] <= 1
        // Therefore this multiplication cannot overflow
        fix16_t product = fix16_mul(dot, u[i]);
        
        fix16_t diff = v[i] - product;
        overflows |= (v[i] ^ product) & (v[i] ^ diff);
        
        v[i] = diff;
    }
    
    if (overflows & SIGNBIT)
    {
        *errors |= FIXMATRIX_OVERFLOW;
    }
}

void mf16_qr_decomposition(mf16 *qt, mf16 *r, const mf16 *matrix, int reorthogonalize)
{
    int i, j, reorth;
    fix16_t dot, norm;
    
    int columns = matrix->columns;
    
    // This uses the modified Gram-Schmidt algorithm.
    // subtract_projection takes advantage of the fact that
    // previous rows have already been normalized.
    
    // We start with q = transpose(matrix), just because it is
    // nicer to work with row vectors in C. Also for mf16_solve
    // it is more useful to get q' anyway.
    if (qt != matrix)
    {
        *qt = *matrix;
    }
    mf16_transpose(qt);
    
    // R is initialized to have square size of cols(A) and zeroed.
    r->columns = columns;
    r->rows = columns;
    r->errors = 0;
    for (j = 0; j < r->rows; j++)
    {
        for (i = 0; i < r->columns; i++)
        {
            r->data[j][i] = 0;
        }
    }
    
    // Now do the actual Gram-Schmidt for the rows.
    for (j = 0; j < qt->rows; j++)
    {
        for (reorth = 0; reorth <= reorthogonalize; reorth++)
        {
            for (i = 0; i < j; i++)
            {
                fix16_t *v = &qt->data[j][0];
                fix16_t *u = &qt->data[i][0];
                
                dot = dotproduct_row(v, u, qt->columns, &qt->errors);
                subtract_projection(v, u, dot, qt->columns, &qt->errors);
                
                r->data[i][j] += dot;
            }
        }
        
        // Normalize the row in q
        dot = dotproduct_row(&qt->data[j][0], &qt->data[j][0], qt->columns, &qt->errors);
        norm = fix16_sqrt(dot);
        r->data[j][j] = norm;
        
        if (norm < 5 && norm > -5)
        {
            // Nearly zero norm, which means that the row
            // was linearly dependent.
            qt->errors |= FIXMATRIX_USEERR;
            continue;
        }
        
        for (i = 0; i < qt->columns; i++)
        {
            // norm >= v[i] for all i, therefore this division
            // doesn't overflow unless norm approaches 0.
            qt->data[j][i] = fix16_div(qt->data[j][i], norm);
        }
    }
    
    r->errors = qt->errors;
}

