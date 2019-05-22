#ifndef _PROJECTOR_H
#define _PROJECTOR_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <mpc.h>
#include <gmp.h>
#include <math.h>
#include <omp.h>
#include <mkl.h>
#include <stdbool.h>
#include <string.h>



// matrix struct
typedef struct {
	int rows;
	MKL_Complex16 *data;
} matrix;

// create (allocate) COMPLEX zero-matrix
matrix* make_matrix(int rows);

// set matrix to zero
void set_zero(matrix *M);

// calculate projector matrix elements
MKL_Complex16 calc_proj_element(double grid_spacing, double xi, double pj, unsigned long int m, unsigned long int n);

// linspace function
void linspace(int xp_lim, int nbr_grid_pts, double grid_spacing, double *xp);

// trace function
long double calc_trace(matrix *M);

// calculate all projectors. Output: array of matrices Pi, and matrix proj_sum
void calc_all_projs(unsigned long int nbr_grid_pts, unsigned long int Nfock, double grid_spacing, double *x, double *p, matrix **Pi, matrix *proj_sum);


void maximum_likelihood(int iterations,long double cutoff, unsigned long int nbr_grid_pts, unsigned long int Nfock,  matrix **Pi, double *histogram, matrix *rho);

matrix *creation_op(int N);

//matrix *displacement(double x, double p, int N);

void displacement(double x, double p, int N, matrix *D);


void calc_projs2(unsigned long int nbr_grid_pts, unsigned long int N, double grid_spacing, double *x, double *p, matrix *rho_thermal, matrix **Pi, matrix *proj_sum);

double calc_thermal(unsigned long int m, unsigned long int nbr_photons);

void test_projs2(unsigned long int nbr_grid_pts, unsigned long int N, double grid_spacing, double *x, double *p, matrix *rho_thermal, matrix *proj_sum);

// when the input state is an ideal thermal state, the projectors will be diagonal
void calc_projs_diag(unsigned long int nbr_grid_pts, unsigned long int N, double grid_spacing, double *x, double *p, matrix *rho_thermal,  MKL_Complex16 **Pi_diag, MKL_Complex16 *Pi_sum);


#endif
