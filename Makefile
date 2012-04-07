# Makefile for running the unittests of libfixmatrix.
CC = gcc

# Basic CFLAGS for debugging
CFLAGS = -g -O0 -Wall -Wextra -Werror -I ../libfixmath -DFIXMATH_NO_CACHE

all: run_unittests

clean:
	rm -f fix16_unittests_???? fixmatrix_unittests


run_unittests: fix16_unittests_noround fix16_unittests_round fixmatrix_unittests
	./fix16_unittests_round > /dev/null
	./fix16_unittests_noround > /dev/null
	./fixmatrix_unittests > /dev/null

fixmatrix_unittests: fixmatrix_unittests.c fixmatrix.c fixmatrix.h ../libfixmath/fix16.c ../libfixmath/fix16_sqrt.c
	$(CC) $(CFLAGS) -o $@ $^

# The fixmatrix unittests are run only in the rounding, overflow detecting
# configuration, but the underflying fix16_base.c is tested also for other
# configurations.
#
# Test naming:
# r = rounding, n = no rounding
# o = overflow detection, n = no overflow detection
# 64 = int64_t math, 32 = int32_t math

run_fix16_unittests: \
	fix16_unittests_ro64 fix16_unittests_no64 \
	fix16_unittests_rn64 fix16_unittests_nn64 \
	fix16_unittests_ro32 fix16_unittests_no32 \
	fix16_unittests_rn32 fix16_unittests_nn32 \
	fix16_unittests_ro8 fix16_unittests_no8 \
	fix16_unittests_rn8 fix16_unittests_nn8
	$(foreach test, $^, \
	echo $(test) && \
	./$(test) > /dev/null && \
	) true

fix16_unittests_no64: DEFINES=-DFIXMATH_NO_ROUNDING
fix16_unittests_rn64: DEFINES=-DFIXMATH_NO_OVERFLOW
fix16_unittests_nn64: DEFINES=-DFIXMATH_NO_ROUNDING -DFIXMATH_NO_OVERFLOW
fix16_unittests_ro32: DEFINES=-DFIXMATH_NO_64BIT
fix16_unittests_no32: DEFINES=-DFIXMATH_NO_ROUNDING -DFIXMATH_NO_64BIT
fix16_unittests_rn32: DEFINES=-DFIXMATH_NO_OVERFLOW -DFIXMATH_NO_64BIT
fix16_unittests_nn32: DEFINES=-DFIXMATH_NO_OVERFLOW -DFIXMATH_NO_ROUNDING -DFIXMATH_NO_64BIT
fix16_unittests_ro8: DEFINES=-DFIXMATH_OPTIMIZE_8BIT
fix16_unittests_no8: DEFINES=-DFIXMATH_NO_ROUNDING -DFIXMATH_OPTIMIZE_8BIT
fix16_unittests_rn8: DEFINES=-DFIXMATH_NO_OVERFLOW -DFIXMATH_OPTIMIZE_8BIT
fix16_unittests_nn8: DEFINES=-DFIXMATH_NO_OVERFLOW -DFIXMATH_NO_ROUNDING -DFIXMATH_OPTIMIZE_8BIT

fix16_unittests_% : fix16_unittests.c fix16_base.c fix16_base.h
	$(CC) $(CFLAGS) $(DEFINES) -o $@ $^ -lm
