# Makefile for running the unittests of libfixmatrix.

# Directory for libfixmath source code
FIXMATH = ../libfixmath

CC = gcc

# Basic CFLAGS for debugging
CFLAGS = -g -O0 -Wall -Wextra -Werror

# Libfixmath rounding is currently broken (Issue 13)
#CFLAGS += -DFIXMATH_NO_ROUNDING
CFLAGS += -DFIXMATH_NO_CACHE

# Include directory for libfixmath
CFLAGS += -I$(FIXMATH)

run_unittests: unittests
	./unittests

unittests: unittests.c fixmatrix.c fixmatrix.h \
	   $(FIXMATH)/fix16.c $(FIXMATH)/fix16_sqrt.c
	$(CC) $(CFLAGS) -o $@ $^

fix16_unittests: fix16_unittests.c fix16_base.c fix16_base.h
	$(CC) $(CFLAGS) -o $@ $^

fix16_unittests_noround: fix16_unittests.c fix16_base.c fix16_base.h
	$(CC) $(CFLAGS) -DFIXMATH_NO_ROUNDING -o $@ $^

fix16_upstream_unittests: fix16_unittests.c $(FIXMATH)/fix16.c fix16_base.h
	$(CC) $(CFLAGS) -o $@ $^
