#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <mkl.h>
#include "functions.h"

//source ../intel/bin/compilervars.sh intel64

int main() {

	unsigned long int i, j, q, nbr_grid_pts, nbr_of_projs, Nfock;
	double alpha = 1;
	double beta = 0;
	double xp_lim_pos, xp_lim_neg, grid_spacing;
	long double trace, pr;

	FILE *f;

	Nfock = 200;

	xp_lim_pos = 22.7207;
	xp_lim_neg = 22.7207;
	grid_spacing = 0.891007;
	nbr_grid_pts = 51;//ceil((xp_lim_pos+xp_lim_neg)/grid_spacing);
	//nbr_grid_pts = 51;
	//grid_spacing = 2.0*xp_lim/nbr_grid_pts;
	nbr_of_projs = nbr_grid_pts*nbr_grid_pts;

	printf("grid spacing %e, grid points %li\n",grid_spacing,nbr_grid_pts);

	// allocate histogram matrix H
	double *H = malloc(nbr_of_projs*sizeof(double));
	// open histogram matrix file
	f = fopen("Bin51_OFF_n1.dat","r");
	// read measurement histogram
	for (i=0; i<nbr_of_projs; i++) {
		fscanf(f, "%le", &H[i]);
	}
	fclose(f);

	//projectors Pi
	matrix **Pi = malloc(nbr_of_projs*sizeof(matrix));

	// allocate phase space grid
	double *x = malloc(nbr_grid_pts*sizeof(double));
	double *p = malloc(nbr_grid_pts*sizeof(double));
	// set phase space grid
	linspace(xp_lim_neg, nbr_grid_pts,grid_spacing, x);
	linspace(xp_lim_neg, nbr_grid_pts,grid_spacing, p);

	// allocate projector sum matrix
	matrix *proj_sum = make_matrix(Nfock,Nfock);

	// calculate all projectors Pi and their sum proj_sum
	calc_all_projs(nbr_grid_pts, Nfock, grid_spacing, x, p, Pi, proj_sum);

	// calculate trace (trace/Nfock should be 1)
	trace = calc_trace(proj_sum);
	printf("trace %Lf\n",trace/Nfock);

	// print projector sum diagonal elements to file (all should be 1)
	f = fopen("diagonals.dat","w");
	for (i=0; i<Nfock; i++)
		fprintf(f,"%f\n",proj_sum->data[i*Nfock + i].real);
	fclose(f);


	// create initial density matrix rho to be (normalized) identity
    matrix *rho = make_matrix(Nfock,Nfock);

	// initialize rho
	for (i=0; i<Nfock; i++) {
		rho->data[i*Nfock + i].real = 1.0/Nfock;
		}
	//rho->data[0].real = 1.0;



	long double cutoff = 1e-20;
	int maxlik_iterations = 50;

	maximum_likelihood(maxlik_iterations, cutoff, nbr_grid_pts, Nfock, Pi, rho);


	f = fopen("reconst_diag.dat","w");
	for (i=0; i<Nfock; i++)
		fprintf(f,"%e\n",rho->data[i*Nfock+i].real);

	fclose(f);
	/********** END MAXLIK *****************/


	// free allocated memory
	for (i=0; i<nbr_of_projs; i++)
		mkl_free(Pi[i]->data);

	mkl_free(Pi);
	free(x);
	free(p);

}
