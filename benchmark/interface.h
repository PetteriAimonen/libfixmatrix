// This file defines the hardware or simulator interface that will be used to
// measure timings and report results.

typedef unsigned long timestamp_t;

// Reset timer/counter/something
void start_timing();

// Return the number of clock cycles passed since start_timing();
timestamp_t end_timing();

// Print the timing values for the given function
void print_timing(const char *function_name, timestamp_t cycles);
