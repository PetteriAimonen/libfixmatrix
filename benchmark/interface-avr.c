#include <avr/io.h>
#include "interface.h"
#include <stdint.h>

void start_timing()
{
    TCCR1B = 1;
    TCNT1 = 0;
}

timestamp_t end_timing()
{
    return TCNT1;
}

#define special_output_port (*((volatile char *)0x20))

void print_timing(const char *function_name, timestamp_t cycles)
{
    uint8_t i = 0;
    while (*function_name)
    {
        special_output_port = *function_name++;
        i++;
    }
    
    while (i++ < 20)
        special_output_port = ' ';
    
    uint8_t first = 1;
    for (i = 0; i < 8; i++)
    {
        uint8_t digit = cycles / 10000000;
        cycles -= digit * 10000000;
        cycles *= 10;
        
        if (!first || digit != 0)
        {
            special_output_port = '0' + digit;
            first = 0;
        } else {
            special_output_port = ' ';
        }
    }
    
    special_output_port = '\n';
}
