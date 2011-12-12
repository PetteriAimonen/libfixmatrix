#include <fix16_base.h>
#include "interface.h"

fix16_t result;
float resultf;
timestamp_t time;

int main()
{
    start_timing();
    time = end_timing();
    print_timing("test", time);
    
    start_timing();
    asm("nop");
    time = end_timing();
    print_timing("test2", time);
    
    
    start_timing();
    result = fix16_div(0x1234102, 0x12352);
    time = end_timing();
    print_timing("fix16_div", time);
    
    start_timing();
    result = fix16_mul(fix16_from_int(123), fix16_from_int(17));
    time = end_timing();
    print_timing("fix16_mul", time);
    print_timing("fix16_mul result", result);
    print_timing("fix16 equ", fix16_mul(0x1FFFF, 0x1FFFF));
    print_timing("test", 262140);
    
    resultf = 3.14159265f;
    start_timing();
    resultf *= 123456.0f;
    time = end_timing();
    print_timing("float mul", time);
    
    start_timing();
    resultf /= 654321.0f;
    time = end_timing();
    print_timing("float div", time);
    
    start_timing();
    result = fix16_sqrt(fix16_from_int(123));
    time = end_timing();
    print_timing("fix16_sqrt", time);
    
    return 0;
}