#ifndef _PROJECTOR_H
#define _PROJECTOR_H
#include <stdlib.h>

int Nfock;

typedef struct {
	int rows;
	int cols;
	double *data;
} matrix;

matrix* make_matrix(int rows, int cols) {
	matrix m;
	m.rows = rows;
	m.cols = cols;
	m.data = calloc(rows*cols, sizeof(double));
}





#endif
