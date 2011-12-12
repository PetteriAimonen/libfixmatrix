#include "interface.h"
#include <stdint.h>
#include <stdio.h>

// This targets an ARM Cortex M3 core using QEmu LM3S6965 emulation.
#define STBASE 0xE000E000
#define STCTRL (*(volatile uint32_t*)(0x010 + STBASE))
#define STRELOAD (*(volatile uint32_t*)(0x014 + STBASE))
#define STCURRENT (*(volatile uint32_t*)(0x018 + STBASE))


void start_timing()
{
     STRELOAD = 0x00FFFFFF;
     STCTRL = 5;
     STCURRENT = 0;
}

timestamp_t end_timing()
{
     return 0x00FFFFFF - STCURRENT - 4;
}

void print_timing(const char *function_name, timestamp_t cycles)
{
    printf("%20s %8d\n", function_name, cycles);
}
