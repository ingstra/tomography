#include "matrix.h"

matrix* make_matrix(int rows, int cols) {
	matrix *m = malloc(sizeof(matrix));
	m->rows = rows;
	m->cols = cols;
	m->data = calloc(rows*cols, sizeof(long double complex));
	return m;
}
