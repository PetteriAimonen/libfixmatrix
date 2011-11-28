#include <stdio.h>
#include "unittests.h"
#include "fixmatrix.h"

void print_matrix(const mf16 *matrix)
{
    if (matrix->errors)
    {
        printf("ERRORS: %d\n", matrix->errors);
    }
    
    int row, column;
    for (row = 0; row < matrix->rows; row++)
    {
        for (column = 0; column < matrix->columns; column++)
        {
            fix16_t value = matrix->data[row][column];
            printf("%9.4f ", fix16_to_float(value));
        }
        printf("\n");
    }
}

fix16_t max_delta(const mf16 *a, const mf16 *b)
{
    fix16_t max = 0;
    int i, j;
    
    if (a->rows != b->rows || a->columns != b->columns ||
        a->errors || b->errors)
    {
        return fix16_max;
    }
    
    for (i = 0; i < a->rows; i++)
    {
        for (j = 0; j < a->columns; j++)
        {
            fix16_t diff = a->data[i][j] - b->data[i][j];
            if (diff < 0) diff = -diff;
            if (diff > max) max = diff;
        }
    }
    
    return max;
}

int main()
{
    int status = 0;
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)},
             {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)},
             {fix16_from_int(7), fix16_from_int(8), fix16_from_int(9)}}};
        mf16 r;
        
        COMMENT("Test 3x3 matrix multiplication");
        mf16_mul(&r, &a, &a);
        TEST(r.errors == 0);
        TEST(r.rows == 3 && r.columns == 3);
        TEST(r.data[0][0] == fix16_from_int(30));
        TEST(r.data[0][1] == fix16_from_int(36));
        TEST(r.data[0][2] == fix16_from_int(42));
        TEST(r.data[1][0] == fix16_from_int(66));
        TEST(r.data[1][1] == fix16_from_int(81));
        TEST(r.data[1][2] == fix16_from_int(96));
        TEST(r.data[2][0] == fix16_from_int(102));
        TEST(r.data[2][1] == fix16_from_int(126));
        TEST(r.data[2][2] == fix16_from_int(150));
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(1000), fix16_from_int(100), fix16_from_int(100)},
             {fix16_from_int(100), fix16_from_int(5), fix16_from_int(6)},
             {fix16_from_int(100), fix16_from_int(8), fix16_from_int(9)}}};
        mf16 r;
        
        COMMENT("Test overflow detection in multiplication");
        
        // Overflow in the multiplication
        mf16_mul(&r, &a, &a);
        TEST(r.errors == FIXMATRIX_OVERFLOW);
        
        // Overflow in summation
        a.data[0][0] = fix16_from_int(150);
        mf16_mul(&r, &a, &a);
        TEST(r.errors == FIXMATRIX_OVERFLOW);
        
        // No overflow
        a.data[0][0] = fix16_from_int(100);
        mf16_mul(&r, &a, &a);
        TEST(r.errors == 0);
    }
    
    {
        mf16 a = {5, 1, 0,
            {{fix16_from_int(101)},
             {fix16_from_int(102)},
             {fix16_from_int(103)},
             {fix16_from_int(104)},
             {fix16_from_int(105)}}};
        mf16 b = {5, 1, 0,
            {{fix16_from_int(51)},
             {fix16_from_int(52)},
             {fix16_from_int(53)},
             {fix16_from_int(54)},
             {fix16_from_int(55)}}};
        
        mf16 atb;
        
        COMMENT("Test transposed multiplication of vectors");
        mf16_mul_t(&atb, &a, &b);
        
        TEST(atb.rows == 1);
        TEST(atb.columns == 1);
        TEST(atb.errors == 0);
        TEST(atb.data[0][0] == fix16_from_int(27305));
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)},
             {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)},
             {fix16_from_int(7), fix16_from_int(8), fix16_from_int(10)}}};
        mf16 q, r, qtq, qr;
        
        COMMENT("Test 3x3 QR-decomposition");
        mf16_qr_decomposition(&q, &r, &a, 3);
        printf("q =\n");
        print_matrix(&q);
        printf("r =\n");
        print_matrix(&r);
        
        mf16_mul_t(&qtq, &q, &q);
        printf("q'q =\n");
        print_matrix(&qtq);
        
        mf16_mul(&qr, &q, &r);
        printf("qr =\n");
        print_matrix(&qr);
        
        fix16_t one = fix16_from_int(1);
        const mf16 identity = {3, 3, 0,
            {{one, 0, 0}, {0, one, 0}, {0, 0, one}}
        };
        TEST(max_delta(&qtq, &identity) < 5);
        TEST(max_delta(&qr, &a) < 5);
    }
    
    {
        mf16 a = {4, 3, 0,
            {{fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)},
             {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)},
             {fix16_from_int(7), fix16_from_int(8), fix16_from_int(9)},
             {fix16_from_int(10), fix16_from_int(11), fix16_from_int(13)},
            }};
        mf16 q, r, qtq, qr;
        
        q = a;
        
        COMMENT("Test 4x3 QR-decomposition with aliasing q = matrix");
        mf16_qr_decomposition(&q, &r, &q, 3);
        printf("q =\n");
        print_matrix(&q);
        printf("r =\n");
        print_matrix(&r);
        
        mf16_mul_t(&qtq, &q, &q);
        printf("q'q =\n");
        print_matrix(&qtq);
        
        mf16_mul(&qr, &q, &r);
        printf("qr =\n");
        print_matrix(&qr);
        
        fix16_t one = fix16_from_int(1);
        const mf16 identity = {3, 3, 0,
            {{one, 0, 0}, {0, one, 0}, {0, 0, one}}
        };
        TEST(max_delta(&qtq, &identity) < 5);
        TEST(max_delta(&qr, &a) < 10);
        
        COMMENT("Test 4x3 QR-decomposition with aliasing r = matrix");
        r = a;
        mf16_qr_decomposition(&q, &r, &r, 3);
        mf16_mul_t(&qtq, &q, &q);
        mf16_mul(&qr, &q, &r);
        TEST(max_delta(&qtq, &identity) < 5);
        TEST(max_delta(&qr, &a) < 10);
    }
    
    {
        mf16 a = {8, 1, 0,
            {{fix16_from_int(1)},
             {fix16_from_int(2)},
             {fix16_from_int(3)},
             {fix16_from_int(4)},
             {fix16_from_int(5)},
             {fix16_from_int(6)},
             {fix16_from_int(7)},
             {fix16_from_int(8)},
            }};
        mf16 q, r, qtq, qr;
        
        COMMENT("Test 8x1 QR-decomposition");
        mf16_qr_decomposition(&q, &r, &a, 3);
        printf("q =\n");
        print_matrix(&q);
        printf("r =\n");
        print_matrix(&r);
        
        mf16_mul_t(&qtq, &q, &q);
        printf("q'q =\n");
        print_matrix(&qtq);
        
        mf16_mul(&qr, &q, &r);
        printf("qr =\n");
        print_matrix(&qr);
        
        const mf16 identity = {1, 1, 0,
            {{fix16_from_int(1)}}
        };
        TEST(max_delta(&qtq, &identity) < 5);
        TEST(max_delta(&qr, &a) < 10);
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)},
             {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)},
             {fix16_from_int(7), fix16_from_int(8), fix16_from_int(10)}}};
        mf16 b = {3, 1, 0,
            {{fix16_from_int(-1)}, {fix16_from_int(-2)}, {fix16_from_int(-3)}}};
        mf16 q, r, x, ax;
        
        COMMENT("Test 3x3 equation solving");
        mf16_qr_decomposition(&q, &r, &a, 3);
        mf16_solve(&x, &q, &r, &b);
        printf("x = \n");
        print_matrix(&x);
        
        mf16_mul(&ax, &a, &x);
        printf("Ax = \n");
        print_matrix(&ax);
        
        TEST(max_delta(&ax, &b) < 10);
    }
    
    {
        mf16 a = {4, 3, 0,
            {{fix16_from_int(31), fix16_from_int(41), fix16_from_int(59)},
             {fix16_from_int(26), fix16_from_int(53), fix16_from_int(58)},
             {fix16_from_int(97), fix16_from_int(93), fix16_from_int(23)},
             {fix16_from_int(84), fix16_from_int(62), fix16_from_int(64)},
            }};
        mf16 b = {4, 1, 0,
            {{fix16_from_int(100)},
             {fix16_from_int(100)},
             {fix16_from_int(100)},
             {fix16_from_int(100)}
            }};
        mf16 q, r, x;
        
        COMMENT("Test 4x3 least squares solving");
        mf16_qr_decomposition(&q, &r, &a, 3);
        mf16_solve(&x, &q, &r, &b);
        printf("x = \n");
        print_matrix(&x);
        
        // Reference result computed using Octave A\b
        mf16 ref = {3, 1, 0,
            {{fix16_from_float(-0.31426f)},
             {fix16_from_float( 1.16055f)},
             {fix16_from_float( 0.90470f)}}};
        TEST(max_delta(&x, &ref) < 5);
    }
    
    if (status != 0)
        fprintf(stdout, "\n\nSome tests FAILED!\n");
    
    return status;

}