#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <mpc.h>
#include <gmp.h>
#include "matrix.h"


long double complex assign_proj_element(double grid_spacing, double xi, double pj, int m, int n) {
	//printf("xi %f xp %f \n",xi,pj);
	//long double complex numerator = cexp(-xi*xi -pj*pj) * cpow(xi + I*pj,m)*cpow(xi - I*pj,n);
	//printf("after xi %f xp %f \n",xi,pj);
	//long double complex denominator = M_PI * csqrt(tgamma(m+1)*tgamma(n+1));
	long double complex integration_approx = grid_spacing*grid_spacing;
	long double complex part1 = cpow(xi + I*pj,m) / csqrt(tgamma(m+1));
	long double complex part2 = cpow(xi - I*pj,n) / csqrt(tgamma(n+1));

	//printf("res %e num %e dom %e %f %f %i %i \n ",creal(integration_approx*numerator/denominator),creal(numerator),creal(denominator),xi,pj,m,n);
	return integration_approx*cexp(-xi*xi -pj*pj) * part1 * part2 / M_PI;

	//return integration_approx*numerator/denominator;
}

// assign gridpoints to array x and p
void linspace(int xp_lim, int nbr_grid_pts, double grid_spacing, double *xp) {
	int i;

	xp[0] = -xp_lim + grid_spacing/2;
	for (i=1; i<nbr_grid_pts; i++) {
		xp[i] = xp[i-1] + grid_spacing;
	}
}

long double calc_trace(matrix *M) {
	int i;
	long double trace = 0;
	int dim = M->rows;

	for (i=0; i < dim; i++) {
		trace += creal(M->data[dim * i + i]);
		//printf("fcn trace %Le %e\n",trace, creal(M->data[dim * i + i]));
	}
	return trace;
}

int main() {

	int i, j, k,nbr_grid_pts, nbr_of_projs, Nfock;
	unsigned long int m, n;
    double xp_lim, grid_spacing;
	long double trace;

	Nfock = 3;

	// projector Pi
	//double *Pi = malloc(Nfock * Nfock * sizeof(double));

	xp_lim = 10.0;
	nbr_grid_pts = 4;
	nbr_of_projs = nbr_grid_pts*nbr_grid_pts;
	grid_spacing = 2.0*xp_lim/nbr_grid_pts;

	printf("grid spacing %e\n",grid_spacing);

	//projectors Pi
	matrix **Pi = malloc(nbr_of_projs*sizeof(matrix));

	// allocate phase space grid
	double *x = malloc(nbr_grid_pts*sizeof(double));
	double *p = malloc(nbr_grid_pts*sizeof(double));
	// set phase space grid
	linspace(xp_lim, nbr_grid_pts,grid_spacing, x);
	linspace(xp_lim, nbr_grid_pts,grid_spacing, p);

	//	for (i=0;i<nbr_grid_pts;i++)
	//	printf("%e\n", x[i]);

	matrix *proj_sum = make_matrix(Nfock,Nfock);

	// allocate 0 matrix
	for (k=0; k<nbr_of_projs; k++)
		Pi[k] = make_matrix(Nfock,Nfock);


	for (i=0; i<nbr_grid_pts; i++) {
		for (j=0; j< nbr_grid_pts; j++) {
			for (m=0; m<Nfock; m++) {
				for (n=0; n<Nfock; n++) {
					//printf("main x %e p %e i %i j %i m %i n %i nbrPi %i element %i \n",x[i],p[j],i,j,m,n,i*nbr_grid_pts + j,m*Nfock + n);
					Pi[i*nbr_grid_pts + j]->data[m*Nfock + n] = assign_proj_element(grid_spacing, x[i], p[j], m, n);
					proj_sum->data[Nfock*m + n] += Pi[i*nbr_grid_pts + j]->data[m*Nfock + n];
				}
			}
		}
	}


	/*

	for (k=0; k<nbr_of_projs; k++) {
		// allocate 0 matrix
		Pi[k] = make_matrix(Nfock,Nfock);
		for (i=0; i<Nfock; i++) {
			for (j=0; j<Nfock; j++) {
				printf("main x %e p %e i %i j %i k %i \n",x[i],p[j],i,j,k);
				Pi[k]->data[Nfock*i + j] = assign_proj_element(grid_spacing, x[i], p[j], i, j);
				proj_sum->data[Nfock*i + j] += Pi[k]->data[Nfock*i + j];
			}
		}
	}

	*/
	m=1;
	trace = calc_trace(proj_sum);
	printf("trace %Lf\n",trace/Nfock);

	FILE *f;
	f = fopen("diagonals.dat","w");
	for (i=0; i<Nfock; i++)
		fprintf(f,"%f\n",creal(proj_sum->data[i*Nfock + i]));


	/*	for (i=0; i<Nfock; i++)
		printf(" %e %e %e \n",creal(Pi[m]->data[Nfock*i + 0]),creal(Pi[m]->data[Nfock*i + 1]),creal(Pi[m]->data[Nfock*i + 2]));
	*/
	printf("elem %e \n", creal(proj_sum->data[0]));


}
