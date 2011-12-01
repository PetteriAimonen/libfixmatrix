# Makefile for running the unittests of libfixmatrix.
CC = gcc

# Basic CFLAGS for debugging
CFLAGS = -g -O0 -Wall -Wextra -Werror

run_unittests: fix16_unittests_noround fix16_unittests_round fixmatrix_unittests
	./fix16_unittests_round > /dev/null
	./fix16_unittests_noround > /dev/null
	./fixmatrix_unittests > /dev/null

clean:
	rm -f fix16_unittests_noround fix16_unittests_round fixmatrix_unittests

fixmatrix_unittests: fixmatrix_unittests.c fixmatrix.c fixmatrix.h fix16_base.c
	$(CC) $(CFLAGS) -o $@ $^

fix16_unittests_round: fix16_unittests.c fix16_base.c fix16_base.h
	$(CC) $(CFLAGS) -o $@ $^ -lm

fix16_unittests_noround: fix16_unittests.c fix16_base.c fix16_base.h
	$(CC) $(CFLAGS) -DFIXMATH_NO_ROUNDING -o $@ $^ -lm
