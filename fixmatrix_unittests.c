#include <stdio.h>
#include "unittests.h"
#include "fixmatrix.h"
#include "fixstring.h"

fix16_t max_delta(const mf16 *a, const mf16 *b)
{
    fix16_t max = 0;
    int i, j;
    
    if (a->rows != b->rows || a->columns != b->columns ||
        a->errors || b->errors)
    {
        return fix16_maximum;
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
        mf16 r, r2;
        
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
        
        COMMENT("Test 3x3 matrix multiplication with aliasing");
        r2 = a;
        mf16_mul(&r2, &r2, &r2);
        TEST(max_delta(&r, &r2) == 0);
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
        
        mf16 atb, abt;
        
        COMMENT("Test mf16_mul_at");
        mf16_mul_at(&atb, &a, &b);
        
        TEST(atb.rows == 1);
        TEST(atb.columns == 1);
        TEST(atb.errors == 0);
        TEST(atb.data[0][0] == fix16_from_int(27305));
        
        COMMENT("Test mf16_mul_bt");
        mf16_mul_bt(&abt, &a, &b);
        
        TEST(abt.rows == 5);
        TEST(abt.columns == 5);
        TEST(abt.errors == 0);
        TEST(abt.data[0][0] == fix16_from_int(101 * 51));
    }
    
    {
        mf16 a = {4, 3, 0,
            {{fix16_from_int(101), fix16_from_int(102), fix16_from_int(103)},
             {fix16_from_int(104), fix16_from_int(105), fix16_from_int(106)},
             {fix16_from_int(107), fix16_from_int(108), fix16_from_int(109)},
             {fix16_from_int(110), fix16_from_int(111), fix16_from_int(112)}
            }};
        mf16 b = {4, 3, 0,
            {{fix16_from_int(-1), fix16_from_int(-2), fix16_from_int(-3)},
             {fix16_from_int(-4), fix16_from_int(-5), fix16_from_int(-6)},
             {fix16_from_int(-7), fix16_from_int(-8), fix16_from_int(-9)},
             {fix16_from_int(-10), fix16_from_int(-11), fix16_from_int(-12)}
            }};
        mf16 ref = {4, 3, 0,
            {{fix16_from_int(100), fix16_from_int(100), fix16_from_int(100)},
             {fix16_from_int(100), fix16_from_int(100), fix16_from_int(100)},
             {fix16_from_int(100), fix16_from_int(100), fix16_from_int(100)},
             {fix16_from_int(100), fix16_from_int(100), fix16_from_int(100)}
            }};
        mf16 r;
        
        COMMENT("Test mf16_add with 4x3 matrices");
        mf16_add(&r, &a, &b);
        TEST(max_delta(&r, &ref) == 0);
        
        COMMENT("Test mf16_add with 4x3 matrices and aliasing");
        r = a;
        mf16_add(&r, &r, &b);
        TEST(max_delta(&r, &ref) == 0);
        
        COMMENT("Test mf16_sub with 4x3 matrices and aliasing");
        mf16_sub(&r, &r, &b);
        TEST(max_delta(&r, &a) == 0);
    }
    
    {
        mf16 a = {5, 1, 0,
            {{fix16_from_int(1)},
             {fix16_from_int(2)},
             {fix16_from_int(20000)},
             {fix16_from_int(-20000)},
             {fix16_from_int(4)}}};
        mf16 b = {5, 1, 0,
            {{fix16_from_int(1)},
             {fix16_from_int(2)},
             {fix16_from_int(20000)},
             {fix16_from_int(3)},
             {fix16_from_int(4)}}};
        mf16 c = {5, 1, 0,
            {{fix16_from_int(1)},
             {fix16_from_int(2)},
             {fix16_from_int(3)},
             {fix16_from_int(20000)},
             {fix16_from_int(4)}}};
        mf16 r;
        
        COMMENT("Test overflow detection in addition");
        mf16_add(&r, &a, &b);
        TEST(r.errors == FIXMATRIX_OVERFLOW);
        mf16_add(&r, &a, &c);
        TEST(r.errors == 0);
        mf16_add(&r, &b, &c);
        TEST(r.errors == 0);
        
        COMMENT("Test overflow detection in subtraction");
        mf16_sub(&r, &a, &c);
        TEST(r.errors == FIXMATRIX_OVERFLOW);
        mf16_sub(&r, &a, &b);
        TEST(r.errors == 0);
        mf16_sub(&r, &b, &c);
        TEST(r.errors == 0);
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)},
             {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)},
             {fix16_from_int(7), fix16_from_int(8), fix16_from_int(9)}}};
        mf16 b = {3, 3, 0,
            {{fix16_from_int(10), fix16_from_int(20), fix16_from_int(30)},
             {fix16_from_int(40), fix16_from_int(50), fix16_from_int(60)},
             {fix16_from_int(70), fix16_from_int(80), fix16_from_int(90)}}};
        mf16 r;
        
        COMMENT("Test multiplication and division by a scalar");
        
        mf16_mul_s(&r, &a, fix16_from_int(10));
        TEST(max_delta(&r, &b) == 0);
        
        mf16_div_s(&r, &r, fix16_from_int(10));
        TEST(max_delta(&r, &a) == 0);
    }
    
    {
        mf16 a = {5, 1, 0,
            {{fix16_from_int(1)},
             {fix16_from_int(2)},
             {fix16_from_int(20000)},
             {fix16_from_int(-20000)},
             {fix16_from_int(4)}}};
        mf16 r;
        
        COMMENT("Test 5x1 transposition");
        mf16_transpose(&r, &a);
        TEST(r.errors == 0);
        TEST(r.rows == 1);
        TEST(r.columns == 5);
        TEST(r.data[0][0] == a.data[0][0]);
        TEST(r.data[0][1] == a.data[1][0]);
        TEST(r.data[0][2] == a.data[2][0]);
        TEST(r.data[0][3] == a.data[3][0]);
        TEST(r.data[0][4] == a.data[4][0]);
        
        COMMENT("Test 1x5 transposition with aliasing");
        mf16_transpose(&r, &r);
        TEST(max_delta(&r, &a) == 0);
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)},
             {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)},
             {fix16_from_int(7), fix16_from_int(8), fix16_from_int(10)}}};
        mf16 q, r, qtq, qr;
        
        COMMENT("Test 3x3 QR-decomposition");
        mf16_qr_decomposition(&q, &r, &a, 1);
        printf("q =\n");
        print_mf16(stdout, &q);
        printf("r =\n");
        print_mf16(stdout, &r);
        
        mf16_mul_at(&qtq, &q, &q);
        printf("q'q =\n");
        print_mf16(stdout, &qtq);
        
        mf16_mul(&qr, &q, &r);
        printf("qr =\n");
        print_mf16(stdout, &qr);
        
        fix16_t one = fix16_from_int(1);
        const mf16 identity = {3, 3, 0,
            {{one, 0, 0}, {0, one, 0}, {0, 0, one}}
        };
        TEST(max_delta(&qtq, &identity) < 5);
        TEST(max_delta(&qr, &a) < 5);
        
        COMMENT("Test 3x3 QR-decomposition without reorthogonalization");
        mf16_qr_decomposition(&q, &r, &a, 0);
        printf("q =\n");
        print_mf16(stdout, &q);
        printf("r =\n");
        print_mf16(stdout, &r);
        
        mf16_mul_at(&qtq, &q, &q);
        printf("q'q =\n");
        print_mf16(stdout, &qtq);
        
        mf16_mul(&qr, &q, &r);
        printf("qr =\n");
        print_mf16(stdout, &qr);
        
        TEST(max_delta(&qtq, &identity) < 10);
        TEST(max_delta(&qr, &a) < 10);
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
        mf16_qr_decomposition(&q, &r, &q, 1);
        printf("q =\n");
        print_mf16(stdout, &q);
        printf("r =\n");
        print_mf16(stdout, &r);
        
        mf16_mul_at(&qtq, &q, &q);
        printf("q'q =\n");
        print_mf16(stdout, &qtq);
        
        mf16_mul(&qr, &q, &r);
        printf("qr =\n");
        print_mf16(stdout, &qr);
        
        fix16_t one = fix16_from_int(1);
        const mf16 identity = {3, 3, 0,
            {{one, 0, 0}, {0, one, 0}, {0, 0, one}}
        };
        TEST(max_delta(&qtq, &identity) < 5);
        TEST(max_delta(&qr, &a) < 10);
        
        COMMENT("Test 4x3 QR-decomposition with aliasing r = matrix");
        r = a;
        mf16_qr_decomposition(&q, &r, &r, 1);
        mf16_mul_at(&qtq, &q, &q);
        mf16_mul(&qr, &q, &r);
        TEST(max_delta(&qtq, &identity) < 5);
        TEST(max_delta(&qr, &a) < 10);
    }
    
    {
        mf16 a = {8, 5, 0,
            {{ 10,     0,     0,     0,   982},
             {  0,     6, -1383,     0,     0},
             {  0,    15,   580,     0,     0},
             {  0,     0,     0,   284,    -3},
             {-56,     0,     0,     0,   284},
             { 64,     0,     0,     0,     0},
             {  0,    64,     0,     0,     0},
             {  0,     0,    64,     0,     0}}
        };
        mf16 q, r, qr;
        
        COMMENT("Test 8x5 QR-decomposition with small values");
        mf16_qr_decomposition(&q, &r, &a, 1);
        
        printf("q =\n");
        print_mf16(stdout, &q);
        printf("r =\n");
        print_mf16(stdout, &r);
        
        TEST(q.errors == 0);
        TEST(r.errors == 0);
        
        mf16_mul(&qr, &q, &r);
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
        mf16_qr_decomposition(&q, &r, &a, 1);
        printf("q =\n");
        print_mf16(stdout, &q);
        printf("r =\n");
        print_mf16(stdout, &r);
        
        mf16_mul_at(&qtq, &q, &q);
        printf("q'q =\n");
        print_mf16(stdout, &qtq);
        
        mf16_mul(&qr, &q, &r);
        printf("qr =\n");
        print_mf16(stdout, &qr);
        
        const mf16 identity = {1, 1, 0,
            {{fix16_from_int(1)}}
        };
        TEST(max_delta(&qtq, &identity) < 10);
        TEST(max_delta(&qr, &a) < 20);
    }
    
    {
        // Carefully chosen to trigger overflow in subtract_projection.
        mf16 a = {3, 3, 0,
            {{fix16_from_int(1), fix16_from_int(32767), fix16_from_int(1)},
             {fix16_from_int(2), fix16_from_int(-32768), fix16_from_int(0)},
             {fix16_from_int(-1), fix16_from_int(32767), fix16_from_int(0)}}};
        mf16 q, r, qtq;
        
        COMMENT("Test intermediate result overflow detection in QR decomp.");
        
        mf16_qr_decomposition(&q, &r, &a, 0);
        mf16_mul_at(&qtq, &q, &q);
        
        printf("q'q =\n");
        print_mf16(stdout, &qtq);
        
        fix16_t one = fix16_from_int(1);
        const mf16 identity = {3, 3, 0,
            {{one, 0, 0}, {0, one, 0}, {0, 0, one}}
        };
        TEST(q.errors == FIXMATRIX_OVERFLOW || max_delta(&qtq, &identity) < 50);
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(535), fix16_from_int(32767), fix16_from_int(1)},
             {fix16_from_int(2), fix16_from_int(23), fix16_from_int(400)},
             {fix16_from_int(324), fix16_from_int(5), fix16_from_int(0)}}};
        mf16 q, r, qtq;
        
        COMMENT("Test large value handling in QR decomp.");
        
        mf16_qr_decomposition(&q, &r, &a, 0);
        mf16_mul_at(&qtq, &q, &q);
        
        printf("q =\n");
        print_mf16(stdout, &q);
        
        printf("q'q =\n");
        print_mf16(stdout, &qtq);
        
        fix16_t one = fix16_from_int(1);
        const mf16 identity = {3, 3, 0,
            {{one, 0, 0}, {0, one, 0}, {0, 0, one}}
        };
        TEST(max_delta(&qtq, &identity) < 50);
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
        mf16_qr_decomposition(&q, &r, &a, 1);
        mf16_solve(&x, &q, &r, &b);
        printf("x = \n");
        print_mf16(stdout, &x);
        
        mf16_mul(&ax, &a, &x);
        printf("Ax = \n");
        print_mf16(stdout, &ax);
        
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
        mf16_qr_decomposition(&q, &r, &a, 1);
        mf16_solve(&x, &q, &r, &b);
        printf("x = \n");
        print_mf16(stdout, &x);
        
        // Reference result computed using Octave A\b
        mf16 ref = {3, 1, 0,
            {{fix16_from_float(-0.31426f)},
             {fix16_from_float( 1.16055f)},
             {fix16_from_float( 0.90470f)}}};
        TEST(max_delta(&x, &ref) < 20);
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(15), fix16_from_int(-12), fix16_from_int(99)},
             {fix16_from_int(42), fix16_from_int(57), fix16_from_int(6)},
             {fix16_from_int(72), fix16_from_int(-8), fix16_from_int(10)}}};
        mf16 b = {3, 2, 0,
            {{fix16_from_int(10), fix16_from_int(-12)},
             {fix16_from_int(20), fix16_from_int(15)},
             {fix16_from_int(30), fix16_from_int(99)}}};
        mf16 q, r, x, ax;
        
        COMMENT("Test 3x3 equation solving with multiple columns");
        mf16_qr_decomposition(&q, &r, &a, 1);
        mf16_solve(&x, &q, &r, &b);
        printf("x = \n");
        print_mf16(stdout, &x);
        
        mf16_mul(&ax, &a, &x);
        printf("Ax = \n");
        print_mf16(stdout, &ax);
        
        // Note: large delta due to large values in matrix.
        // This is one of the shortcomings of fixed point format.
        TEST(max_delta(&ax, &b) < 100);
    }
    
    {
        mf16 a = {4, 4, 0,
            {{fix16_from_int(7), fix16_from_int(-11), fix16_from_int(80), fix16_from_int(15)},
             {fix16_from_int(11), fix16_from_int(-59), fix16_from_int(57), fix16_from_int(72)},
             {fix16_from_int(79), fix16_from_int(57), fix16_from_int(-8), fix16_from_int(24)},
             {fix16_from_int(-23), fix16_from_int(32), fix16_from_int(0), fix16_from_int(56)},
            }};
        fix16_t one = fix16_from_int(1);
        mf16 identity = {4, 4, 0,
            {{one,0,0,0}, {0,one,0,0}, {0,0,one,0}, {0,0,0,one}}};
        mf16 q, r, result, inv_a;
        
        COMMENT("Test 4x4 matrix inversion");
        mf16_qr_decomposition(&q, &r, &a, 1);
        
        printf("q =\n");
        print_mf16(stdout, &q);
        printf("r =\n");
        print_mf16(stdout, &r);
        
        mf16_solve(&inv_a, &q, &r, &identity);
        printf("inv(a) =\n");
        print_mf16(stdout, &inv_a);
        
        mf16_mul(&result, &a, &inv_a);
        printf("a*inv(a) =\n");
        print_mf16(stdout, &result);
        TEST(max_delta(&result, &identity) < 100);
        
        mf16_mul(&result, &inv_a, &a);
        printf("inv(a)*a =\n");
        print_mf16(stdout, &result);
        TEST(max_delta(&result, &identity) < 100);
        
        COMMENT("Test 4x4 matrix inversion with aliasing dest = matrix");
        inv_a = identity;
        mf16_solve(&inv_a, &q, &r, &inv_a);
        printf("inv(a) =\n");
        print_mf16(stdout, &inv_a);
        mf16_mul(&result, &a, &inv_a);
        TEST(max_delta(&result, &identity) < 100);
        
        COMMENT("Test 4x4 matrix inversion with aliasing dest = q");
        inv_a = q;
        mf16_solve(&inv_a, &inv_a, &r, &identity);
        printf("inv(a) =\n");
        print_mf16(stdout, &inv_a);
        mf16_mul(&result, &a, &inv_a);
        TEST(max_delta(&result, &identity) < 100);
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(66), fix16_from_int(78), fix16_from_int(90)},
             {fix16_from_int(78), fix16_from_int(93), fix16_from_int(108)},
             {fix16_from_int(90), fix16_from_int(108), fix16_from_int(126)}}};
        mf16 l, l2, llt;
        
        COMMENT("Test 3x3 Cholesky factorization");
    
        mf16_cholesky(&l, &a);
        mf16_mul_bt(&llt, &l, &l);
        
        printf("L =\n");
        print_mf16(stdout, &l);
        
        printf("L L' =\n");
        print_mf16(stdout, &llt);
        
        TEST(max_delta(&a, &llt) < 10);
        
        COMMENT("Test Cholesky factorization with aliasing");
        l2 = a;
        mf16_cholesky(&l2, &l2);
        TEST(max_delta(&l, &l2) == 0);
    }
    
    {
        mf16 a = {3, 3, 0,
            {{fix16_from_int(2), fix16_from_int(-1), fix16_from_int(0)},
             {fix16_from_int(-1), fix16_from_int(2), fix16_from_int(-1)},
             {fix16_from_int(0), fix16_from_int(-1), fix16_from_int(2)}}};
        mf16 l;
             
        fix16_t one = fix16_from_int(1);
        mf16 identity = {3, 3, 0,
            {{one,0,0}, {0,one,0}, {0,0,one}}};
        
        COMMENT("Test mf16_invert_lt");
        
        printf("A =\n");
        print_mf16(stdout, &a);
        
        mf16_cholesky(&l, &a);
        printf("L =\n");
        print_mf16(stdout, &l);
        
        mf16_invert_lt(&l, &l);
        
        printf("inv A =\n");
        print_mf16(stdout, &l);
        
        mf16_mul(&a, &l, &a);
        printf("inv A * A =\n");
        print_mf16(stdout, &a);
        
        TEST(max_delta(&a, &identity) < 10);
    }
        
    if (status != 0)
        fprintf(stdout, "\n\nSome tests FAILED!\n");
    
    return status;

}
