CC = gcc
CFLAGS = -O2 -Wall -S
BLASFLAGS = -L/usr/lib64/atlas -lsatlas

all: matrixadd matrixmul matrixmul_blas

matrixadd: matrixadd.c matrix.h
	$(CC) $(CFLAGS) -o matrixadd matrixadd.c

matrixmul: matrixmul.c matrix.h
	$(CC) $(CFLAGS) -o matrixmul matrixmul.c

matrixmul_blas: matrixmul_blas.c matrix.h
	$(CC) $(CFLAGS) -o matrixmul_blas matrixmul_blas.c $(BLASFLAGS)

clean:
	rm -f matrixadd matrixmul matrixmul_blas
