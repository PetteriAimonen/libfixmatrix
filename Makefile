# Makefile for running the unittests of libfixmatrix.
CC = gcc

# Basic CFLAGS for debugging
CFLAGS = -g -O0 -Wall -Wextra -Werror -I libfixmath -DFIXMATH_NO_CACHE

COMMON = fixarray.c fixstring.c libfixmath/fix16.c libfixmath/fix16_sqrt.c

all: run_unittests

clean:
	rm -f fixmatrix_unittests

run_unittests: fixmatrix_unittests fixmatrix_unittests_32bit fixvector3d_unittests
	./fixmatrix_unittests > /dev/null
	./fixmatrix_unittests_32bit > /dev/null
	./fixvector3d_unittests > /dev/null

fixmatrix_unittests: fixmatrix_unittests.c fixmatrix.c fixmatrix.h $(COMMON)
	$(CC) $(CFLAGS) -o $@ $^

fixmatrix_unittests_32bit: fixmatrix_unittests.c fixmatrix.c fixmatrix.h $(COMMON)
	$(CC) $(CFLAGS) -DFIXMATH_NO_64BIT -o $@ $^

fixvector3d_unittests: fixvector3d_unittests.c fixvector3d.c fixvector3d.h $(COMMON)
	$(CC) $(CFLAGS) -o $@ $^

libfixmath/%:
	@echo "Downloading a copy of libfixmath..."
	svn co http://libfixmath.googlecode.com/svn/trunk/libfixmath

