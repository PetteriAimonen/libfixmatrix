# Makefile for running the unittests of libfixmatrix.
CC = gcc

# Basic CFLAGS for debugging
CFLAGS = -g -O0 -Wall -Wextra -Werror -I libfixmath -DFIXMATH_NO_CACHE

all: run_unittests

clean:
	rm -f fixmatrix_unittests

run_unittests: fixmatrix_unittests
	./fixmatrix_unittests > /dev/null

fixmatrix_unittests: fixmatrix_unittests.c fixmatrix.c fixmatrix.h libfixmath/fix16.c libfixmath/fix16_sqrt.c
	$(CC) $(CFLAGS) -o $@ $^

libfixmath/%:
	@echo "Downloading a copy of libfixmath..."
	svn co http://libfixmath.googlecode.com/svn/trunk/libfixmath

