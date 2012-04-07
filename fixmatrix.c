#include "fixmatrix.h"

/****************************
 * Initialization functions *
 ****************************/

void mf16_fill(mf16 *dest, fix16_t value)
{
    int row, column;
    dest->errors = 0;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = value;
        }
    }
}

void mf16_fill_diagonal(mf16 *dest, fix16_t value)
{
    int row;
    
    mf16_fill(dest, 0);
    
    for (row = 0; row < dest->rows; row++)
    {
        dest->data[row][row] = value;
    }
}


/*********************************
 * Operations between 2 matrices *
 *********************************/

// Calculates the dotproduct of two vectors of size n.
// If overflow happens, sets flag in errors
fix16_t dotproduct(const fix16_t *a, uint8_t a_stride, const fix16_t *b, uint8_t b_stride, uint8_t n, uint8_t *errors)
{
    fix16_t sum = 0;
    
    while (n--)
    {
        // Compute result
        fix16_t product = fix16_mul(*a, *b);
        sum = fix16_add(sum, product);
        
        if (sum == fix16_overflow || product == fix16_overflow)
            *errors |= FIXMATRIX_OVERFLOW;
        
        // Go to next item
        a += a_stride;
        b += b_stride;
    }
    
    return sum;
}

void mf16_mul(mf16 *dest, const mf16 *a, const mf16 *b)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp destination.
    mf16 tmp;
    mf16 *realdest = dest;
    if (dest == a || dest == b)
        dest = &tmp;
    
    dest->errors = a->errors | b->errors;
    
    if (a->columns != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
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
    
    if (dest != realdest)
        *realdest = *dest;
}

// Multiply transpose of at with b
void mf16_mul_at(mf16 *dest, const mf16 *at, const mf16 *b)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp destination.
    mf16 tmp;
    mf16 *realdest = dest;
    if (dest == at || dest == b)
        dest = &tmp;
    
    dest->errors = at->errors | b->errors;
    
    if (at->rows != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
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
    
    if (dest != realdest)
        *realdest = *dest;
}

void mf16_mul_bt(mf16 *dest, const mf16 *a, const mf16 *bt)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp destination.
    mf16 tmp;
    mf16 *realdest = dest;
    if (dest == a || dest == bt)
        dest = &tmp;
    
    dest->errors = a->errors | bt->errors;
    
    if (a->columns != bt->columns)
        dest->errors |= FIXMATRIX_DIMERR;
    
    dest->rows = a->rows;
    dest->columns = bt->rows;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = dotproduct(
                &a->data[row][0], 1,
                &bt->data[column][0], 1,
                a->columns, &dest->errors);
        }
    }
    
    if (dest != realdest)
        *realdest = *dest;
}

static void mf16_addsub(mf16 *dest, const mf16 *a, const mf16 *b, uint8_t add)
{
    int row, column;
    
    dest->errors = a->errors | b->errors;
    if (a->columns != b->columns || a->rows != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
    dest->rows = a->rows;
    dest->columns = a->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            fix16_t sum;
            if (add)
                sum = fix16_add(a->data[row][column], b->data[row][column]);
            else
                sum = fix16_sub(a->data[row][column], b->data[row][column]);
                
            if (sum == fix16_overflow)
                dest->errors |= FIXMATRIX_OVERFLOW;
            
            dest->data[row][column] = sum;
        }
    }
}

void mf16_add(mf16 *dest, const mf16 *a, const mf16 *b)
{
    mf16_addsub(dest, a, b, 1);
}

void mf16_sub(mf16 *dest, const mf16 *a, const mf16 *b)
{
    mf16_addsub(dest, a, b, 0);
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
    
    dest->rows = matrix->rows;
    dest->columns = matrix->columns;
    dest->errors = matrix->errors;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            fix16_t value = matrix->data[row][column];
            fix16_t product = fix16_mul(value, scalar);
            
            if (product == fix16_overflow)
                dest->errors |= FIXMATRIX_OVERFLOW;
            
            dest->data[row][column] = product;
        }
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
    while (n--)
    {
        // For unit vector u, u[i] <= 1
        // Therefore this multiplication cannot overflow
        fix16_t product = fix16_mul(dot, *u);
        
        // Overflow here is rare, but possible.
        fix16_t diff = fix16_sub(*v, product);
        
        if (diff == fix16_overflow)
            *errors |= FIXMATRIX_OVERFLOW;
        
        *v = diff;
        
        v += FIXMATRIX_MAX_SIZE;
        u += FIXMATRIX_MAX_SIZE;
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
        uint8_t dp_errors = 0;
        dot = dotproduct(&q->data[0][j], stride, &q->data[0][j], stride,
                         n, &dp_errors);
        norm = fix16_sqrt(dot);
        r->data[j][j] = norm;
        
        // Dot product may overflow if the values are larger than 256.
        // If this happens, prescale them before calculating dot product.
        if (dp_errors)
        {
            fix16_t norm_scaler = fix16_from_int(256);
            for (i = 0; i < n; i++) {
                q->data[i][j] = fix16_div(q->data[i][j], norm_scaler);
            }
            
            dot = dotproduct(&q->data[0][j], stride, &q->data[0][j], stride,
                         n, &q->errors);
            
            norm = fix16_sqrt(dot);
            r->data[j][j] = fix16_mul(norm, norm_scaler);
        }
        
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
    
    if (r->columns != r->rows || r->columns != q->columns || r == dest)
    {
        dest->errors |= FIXMATRIX_USEERR;
        return;
    }
    
    // Ax=b <=> QRx=b <=> Q'QRx=Q'b <=> Rx=Q'b
    // Q'b is calculated directly and x is then solved row-by-row.
    mf16_mul_at(dest, q, matrix);
    
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
                value = fix16_sub(value, product);
                
                if (product == fix16_overflow || value == fix16_overflow)
                {
                    dest->errors |= FIXMATRIX_OVERFLOW;
                }
            }
            
            // Now value = R_ij x_i <=> x_i = value / R_ij
            fix16_t divider = r->data[row][row];
            if (divider == 0)
            {
                dest->errors |= FIXMATRIX_SINGULAR;
                dest->data[row][column] = 0;
                continue;
            }
            
            fix16_t result = fix16_div(value, divider);
            dest->data[row][column] = result;
            
            if (result == fix16_overflow)
            {
                dest->errors |= FIXMATRIX_OVERFLOW;
            }
        }
    }
}


