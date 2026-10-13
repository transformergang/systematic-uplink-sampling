# Native build of the portable integer decision code and its bit-exactness tests.
# Needs only a C99 compiler:  make test   (CC=cc by default; override with CC=gcc, CC=clang, CC="zig cc", ...)
CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra
INC      = -Ialg -Ivectors
ALG      = $(wildcard alg/*.c)

all: test

build/test_native_madow: $(ALG) tests/test_native_madow_main.c $(wildcard alg/*.h) vectors/vectors.h vectors/vectors_madow.h
	mkdir -p build
	$(CC) $(CFLAGS) $(INC) -o $@ $(ALG) tests/test_native_madow_main.c

# 34 vectors of the closed-loop family (independent, EXCL, EXCL with online block starts, duty cycle, trigger)
# + 32 vectors of the open-loop single-start rule + 4 self-checks; the last line must start with PASS.
test: build/test_native_madow
	./build/test_native_madow | tee build/test_native_madow.log | tail -n 8
	grep -q '^PASS' build/test_native_madow.log

# compile-only check of the band-level score
score-check: score/score.c score/score.h score/tables.h
	mkdir -p build
	$(CC) $(CFLAGS) -Iscore -c score/score.c -o build/score.o

clean:
	rm -rf build
.PHONY: all test score-check clean
