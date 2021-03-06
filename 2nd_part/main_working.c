#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <mkl.h>
#include "functions.h"

//source ../intel/bin/compilervars.sh intel64

int main() {

	unsigned long int i, j, q, nbr_grid_pts, nbr_of_projs, Nfock;
	double xp_lim_pos, xp_lim_neg, grid_spacing;
	long double trace;

	FILE *f;

	Nfock = 300;

	//Bin51_OFF_n1.dat  Rows: 51 Start: -22.7207 Delta: 0.891007
	//
	xp_lim_pos = 22.156;
	xp_lim_neg = xp_lim_pos;
	grid_spacing = 0.312056;
	//nbr_grid_pts = ceil((xp_lim_pos+xp_lim_neg)/grid_spacing);
	nbr_grid_pts = 143;
	//grid_spacing = 2.0*xp_lim/nbr_grid_pts;
	nbr_of_projs = nbr_grid_pts*nbr_grid_pts;

	printf("grid spacing %e, grid points %li\n",grid_spacing,nbr_grid_pts);

	// allocate histogram matrix H
	double *H = malloc(nbr_of_projs*sizeof(double));
	// open histogram matrix file
	f = fopen("Bin143_start_-22.156_delta_0.312056_OFF.txt","r");
	// read measurement histogram
	for (i=0; i<nbr_of_projs; i++) {
		fscanf(f, "%le", &H[i]);
	}
	fclose(f);

	/*	for (i=0; i<nbr_of_projs; i++)
		printf("hist %le \n",H[i]);
	*/
	//projectors Pi
	matrix **Pi = malloc(nbr_of_projs*sizeof(matrix));

	// allocate phase space grid
	double *x = malloc(nbr_grid_pts*sizeof(double));
	double *p = malloc(nbr_grid_pts*sizeof(double));
	// set phase space grid
	linspace(xp_lim_neg, nbr_grid_pts,grid_spacing, x);
	linspace(xp_lim_neg, nbr_grid_pts,grid_spacing, p);

	// allocate first projector sum matrix
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


	// create initial density matrix rho
    matrix *rho = make_matrix(Nfock,Nfock);
	// initialize rho
	for (i=0; i<Nfock; i++) {
		rho->data[i*Nfock + i].real = 1.0/Nfock;
		}

	// maxlik
	long double cutoff = 1e-20;
	int maxlik_iterations = 100;
	maximum_likelihood(maxlik_iterations, cutoff, nbr_grid_pts, Nfock, Pi, H, rho);

	// print result
	f = fopen("reconst_diag.dat","w");
	for (i=0; i<Nfock; i++)
		fprintf(f,"%e\n",rho->data[i*Nfock+i].real);
	fclose(f);


	// free allocated memory
	for (i=0; i<nbr_of_projs; i++)
		mkl_free(Pi[i]->data);

	mkl_free(Pi);
	free(x);
	free(p);
	free(H);

}
