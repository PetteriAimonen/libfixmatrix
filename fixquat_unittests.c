#include <stdio.h>
#include "unittests.h"
#include "fixquat.h"
#include "fixstring.h"

fix16_t max_delta(const qf16 *a, const qf16 *b)
{
    fix16_t max = 0;
    max = fix16_max(max, fix16_abs(a->a - b->a));
    max = fix16_max(max, fix16_abs(a->b - b->b));
    max = fix16_max(max, fix16_abs(a->c - b->c));
    max = fix16_max(max, fix16_abs(a->d - b->d));
    return max;
}

int main()
{
    int status = 0;
    
    {
        qf16 a = {fix16_from_int(1), fix16_from_int(2), fix16_from_int(3), fix16_from_int(4)};
        qf16 b = {fix16_from_int(5), fix16_from_int(6), fix16_from_int(7), fix16_from_int(8)};
        qf16 c;
        
        COMMENT("Test basic arithmetic");
        
        qf16_conj(&c, &a);
        TEST(c.a == fix16_from_int(1));
        TEST(c.b == fix16_from_int(-2));
        TEST(c.c == fix16_from_int(-3));
        TEST(c.d == fix16_from_int(-4));
        
        qf16_mul(&c, &a, &b);
        TEST(c.a == fix16_from_int(-60));
        TEST(c.b == fix16_from_int(12));
        TEST(c.c == fix16_from_int(30));
        TEST(c.d == fix16_from_int(24));
        
        qf16_mul_s(&c, &a, fix16_from_int(2));
        TEST(c.a == fix16_from_int(2));
        TEST(c.b == fix16_from_int(4));
        TEST(c.c == fix16_from_int(6));
        TEST(c.d == fix16_from_int(8));
        
        qf16_div_s(&c, &c, fix16_from_int(2));
        TEST(max_delta(&a, &c) == 0);
        
        COMMENT("Test norm");
        TEST(fix16_abs(qf16_norm(&a) - fix16_from_float(5.477225f)) < 2);
        
        COMMENT("Test normalize");
        qf16_normalize(&c, &a);
        TEST(fix16_abs(c.a - fix16_from_float(0.182574f)) < 2);
        TEST(fix16_abs(c.b - fix16_from_float(0.365148f)) < 2);
        TEST(fix16_abs(c.c - fix16_from_float(0.547722f)) < 2);
        TEST(fix16_abs(c.d - fix16_from_float(0.730297f)) < 2);
        
    }
    
    {
        qf16 rot = {fix16_from_float(0.7071), fix16_from_float(0.7071), 0, 0};
        mf16 matrix;
        
        qf16_to_matrix(&matrix, &rot);
        
        print_qf16(stdout, &rot);
        printf("\n\n");
        print_mf16(stdout, &matrix);
        
        TEST(matrix.data[0][0] == fix16_from_int(1));
        TEST(matrix.data[1][0] == fix16_from_int(0));
        TEST(matrix.data[2][0] == fix16_from_int(0));
        TEST(matrix.data[0][1] == fix16_from_int(0));
        TEST(matrix.data[1][1] == fix16_from_int(0));
        TEST(matrix.data[2][1] == fix16_from_int(1));
        TEST(matrix.data[0][2] == fix16_from_int(0));
        TEST(matrix.data[1][2] == fix16_from_int(-1));
        TEST(matrix.data[2][2] == fix16_from_int(0));
    }
    
    if (status != 0)
        fprintf(stdout, "\n\nSome tests FAILED!\n");
    
    return status;
}