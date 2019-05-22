#ifndef _MATRIX_EXPONENTIAL_H
#define _MATRIX_EXPONENTIAL_H
#include "functions.h"

/*
  Complex functions.
*/
void expm( matrix *a, matrix *e);
double complex *c8mat_expm1 ( int n, double complex a[] );
/*
  Real functions.
*/
double *r8mat_expm1 ( int n, double a[] );
double *r8mat_expm2 ( int n, double a[] );
double *r8mat_expm3 ( int n, double a[] );

#endif
