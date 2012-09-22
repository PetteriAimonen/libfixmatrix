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

fix16_t max_delta_abs(const qf16 *a, const qf16 *b)
{
    qf16 a_inv = {-a->a, -a->b, -a->c, -a->d};
    return fix16_min(max_delta(a, b), max_delta(&a_inv, b));
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
        COMMENT("Test quaternion power");
        qf16 a, b, a_sq, result;
        v3d axis = {F16(0.8), F16(0.6), 0};
        
        qf16_from_axis_angle(&a, &axis, F16(1.0));
        qf16_from_axis_angle(&b, &axis, F16(1.6));
        
        qf16_mul(&a_sq, &a, &a);
        qf16_pow(&result, &a, F16(2.0));
        printf("(");
        print_qf16(stdout, &a);
        printf(") ^ 2 = ");
        print_qf16(stdout, &result);
        printf("   should be   ");
        print_qf16(stdout, &a_sq);
        printf("\n");
        TEST(max_delta(&result, &a_sq) < F16(0.01));
        
        qf16_pow(&result, &a, F16(1.6));
        printf("(");
        print_qf16(stdout, &a);
        printf(") ^ 1.6 = ");
        print_qf16(stdout, &result);
        printf("\n");
        TEST(max_delta(&result, &b) < F16(0.01));
        
        qf16_pow(&result, &b, F16(0.625));
        printf("(");
        print_qf16(stdout, &b);
        printf(") ^ 0.625 = ");
        print_qf16(stdout, &result);
        printf("\n");
        TEST(max_delta(&result, &a) < F16(0.01));
    }
    
    {
        COMMENT("Test qf16 averaging");
        qf16 a = {F16(1), 0, 0, 0}, b = {-F16(1), 0, 0, 0};
        qf16 result, expected;
        
        qf16_avg(&result, &a, &b, F16(0.2));
        printf("avg = ");
        print_qf16(stdout, &result);
        printf("\n");
        TEST(max_delta_abs(&result, &a) < 2);
        
        a.c = F16(-3);
        b.c = F16(3);
        qf16_normalize(&a, &a);
        qf16_normalize(&b, &b);
        qf16_avg(&result, &a, &b, F16(0.9));
        printf("avg = ");
        print_qf16(stdout, &result);
        printf("\n");
        TEST(max_delta_abs(&result, &a) < 2);
        
        v3d axis = {0, 0, F16(1)};
        qf16_from_axis_angle(&a, &axis, F16(1.0));
        qf16_from_axis_angle(&b, &axis, F16(2.0));
        qf16_from_axis_angle(&expected, &axis, F16(1.5));
        qf16_avg(&result, &a, &b, F16(0.5));
        print_qf16(stdout, &result);
        printf(" vs. ");
        print_qf16(stdout, &expected);
        printf("\n");
        TEST(max_delta_abs(&result, &expected) < 2);
    }
    
    {
        COMMENT("Test qf16 -> rotation matrix");
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
    
    {
        COMMENT("Test qf16_rotate");
        qf16 rot = {F16(0.5), F16(0.5), F16(0.5), F16(0.5)};
        v3d input = {F16(1), F16(2), F16(3)};
        v3d output;
        
        qf16_rotate(&output, &rot, &input);
        TEST(output.x == F16(3));
        TEST(output.y == F16(1));
        TEST(output.z == F16(2));
    }
    
    if (status != 0)
        fprintf(stdout, "\n\nSome tests FAILED!\n");
    
    return status;
}