#include <stdio.h>
#include "unittests.h"
#include "fixvector3d.h"

fix16_t max_delta(const v3d *a, const v3d *b)
{
    fix16_t max = 0;
    max = fix16_max(max, fix16_abs(a->x - b->x));
    max = fix16_max(max, fix16_abs(a->y - b->y));
    max = fix16_max(max, fix16_abs(a->z - b->z));
    return max;
}

int main()
{
    int status = 0;
    
    {
        v3d a = {fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)};
        v3d b = {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)};
        v3d c;
        
        COMMENT("Test basic arithmetic");
        v3d_add(&c, &a, &b);
        TEST(c.x == fix16_from_int(5));
        TEST(c.y == fix16_from_int(7));
        TEST(c.z == fix16_from_int(9));
        
        v3d_sub(&c, &a, &b);
        TEST(c.x == fix16_from_int(-3));
        TEST(c.y == fix16_from_int(-3));
        TEST(c.z == fix16_from_int(-3));
        
        v3d_mul_s(&c, &a, fix16_from_int(2));
        TEST(c.x == fix16_from_int(2));
        TEST(c.y == fix16_from_int(4));
        TEST(c.z == fix16_from_int(6));
        
        v3d_div_s(&c, &c, fix16_from_int(2));
        TEST(max_delta(&c, &a) == 0);
    }
    
    {
        v3d medium = {fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)};
        v3d small = {10, 20, 30};
        v3d large = {fix16_from_int(10000), fix16_from_int(15000), fix16_from_int(20000)};
        
        COMMENT("Test v3d_norm");
        TEST(fix16_abs(v3d_norm(&medium) - fix16_from_float(3.741657f)) < 2);
        TEST(fix16_abs(v3d_norm(&small) - fix16_from_float(0.0005709f)) < 2);
        TEST(fix16_abs(v3d_norm(&large) - fix16_from_float(26925.824f)) < 2);
    }
    
    {
        v3d a = {fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)};
        v3d b = {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)};
        
        COMMENT("Test v3d_dot");
        TEST(v3d_dot(&a, &b) == fix16_from_int(32));
    }
    
    {
        v3d a = {fix16_from_int(3), fix16_from_int(-3), fix16_from_int(1)};
        v3d b = {fix16_from_int(4), fix16_from_int(9), fix16_from_int(2)};
        v3d c;
        v3d expected = {fix16_from_int(-15), fix16_from_int(-2), fix16_from_int(39)};
        
        COMMENT("Test v3d_cross");
        v3d_cross(&c, &a, &b);
        TEST(max_delta(&expected, &c) < 2);
    }
    
    if (status != 0)
        fprintf(stdout, "\n\nSome tests FAILED!\n");
    
    return status;
}