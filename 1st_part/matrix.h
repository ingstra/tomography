#ifndef _PROJECTOR_H
#define _PROJECTOR_H
#include <stdlib.h>
#include <complex.h>

int Nfock;

typedef struct {
	int rows;
	int cols;
	double complex *data;
} matrix;

matrix* make_matrix(int rows, int cols);





#endif
