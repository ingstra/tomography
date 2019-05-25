CC = gcc
LIBS = -lm -lmpc -lgmp -lblas -llapack -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp
CFLAGS = -O3 -Wall -g -fopenmp -march=native -ffast-math
LDFLAGS = -L/home/ingrid/intel/compilers_and_libraries_2019.2.187/linux/mkl/lib/intel64
OBJS = main.o functions.o matrix_exponential.o c8lib.o r8lib.o


all: main

main: $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

main.o: main.c functions.h

functions.o: functions.h

matrix_exponential.o: matrix_exponential.c

c8lib.o: c8lib.c
r8lib.o: c8lib.c

clean:
	rm -f main $(OBJS)
