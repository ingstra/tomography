#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <mkl.h>
#include "functions.h"
#include "matrix_exponential.h"


//source ~/intel/bin/compilervars.sh intel64

int main() {
	unsigned long int i, j, q, nbr_grid_pts, nbr_of_projs, Nfock, nbr_photons;
	double xp_lim_pos, xp_lim_neg, grid_spacing;
	long double trace;

	FILE *f;

	Nfock = 500;

	mkl_set_dynamic(0);
	mkl_set_num_threads(28);

	//Bin51_OFF_n1.dat  Rows: 51 Start: -22.7207 Delta: 0.891007
	//
	xp_lim_pos = 22.156;
	xp_lim_neg = xp_lim_pos;
	grid_spacing =0.312056; //0.1196;
	//nbr_grid_pts = ceil((xp_lim_pos+xp_lim_neg)/grid_spacing);
	nbr_grid_pts = 143;
	//grid_spacing = (xp_lim_pos+xp_lim_neg)/nbr_grid_pts;
	nbr_of_projs = nbr_grid_pts*nbr_grid_pts;

	printf("grid spacing %e, grid points %li\n",grid_spacing,nbr_grid_pts);

	// allocate histogram matrix H
	double *H = malloc(nbr_of_projs*sizeof(double));
	// open histogram matrix file
	f = fopen("Bin143_start_-22.156_delta_0.312056_ON.txt","r");
	// read measurement histogram
	for (i=0; i<nbr_of_projs; i++) {
		fscanf(f, "%le", &H[i]);
	}
	fclose(f);

	/*	for (i=0; i<nbr_of_projs; i++)
		printf("hist %le \n",H[i]);
	*/
	//projectors Pi
	//matrix **Pi = malloc(nbr_of_projs*sizeof(matrix));

	// allocate phase space grid
	double *x = malloc(nbr_grid_pts*sizeof(double));
	double *p = malloc(nbr_grid_pts*sizeof(double));
	// set phase space grid
	linspace(xp_lim_neg, nbr_grid_pts,grid_spacing, x);
	linspace(xp_lim_neg, nbr_grid_pts,grid_spacing, p);

	// allocate first projector sum matrix
	//	matrix *proj_sum = make_matrix(Nfock,Nfock);
	// calculate all projectors Pi and their sum proj_sum
	//calc_all_projs(nbr_grid_pts, Nfock, grid_spacing, x, p, Pi, proj_sum);
	// calculate trace (trace/Nfock should be 1)
	//trace = calc_trace(proj_sum);
	// 	printf("trace %Lf\n",trace/Nfock);
	// print projector sum diagonal elements to file (all should be 1)
	/*	f = fopen("diagonals.dat","w");
		for (i=0; i<Nfock; i++)
		fprintf(f,"%f\n",proj_sum->data[i*Nfock + i].real);
		fclose(f);
	*/



	// maxlik
	long double cutoff = 1e-20;
	int maxlik_iterations = 300;
	//maximum_likelihood(maxlik_iterations, cutoff, nbr_grid_pts, Nfock, Pi, H, rho);

	// print diagonal result
	/*	f = fopen("reconst_diag.dat","w");
		for (i=0; i<Nfock; i++)
		fprintf(f,"%e\n",rho->data[i*Nfock+i].real);
		fclose(f);
	*/
	// print entire state
	/*	f = fopen("reconst_thermal_state.dat", "w");
		fprintf(f,"Nfock %li\n",Nfock);
		for (i=0; i<Nfock*Nfock; i++)
		fprintf(f,"%f %f\n",rho->data[i].real,rho->data[i].imag);
		fclose(f);*/


	/************ second part ********************/

	/*			matrix *D = make_matrix(Nfock,Nfock);
				displacement(0.25,0,Nfock,D);
				for (i=0;i<3;i++)
				printf("%f %f %f \n",D->data[i*3].real,D->data[i*3+1].real,D->data[i*3+2].real);*/


	nbr_photons = 35;
	// create density matrix rho
    matrix *rho = make_matrix(Nfock);
	// set rho to a thermal state
	set_zero(rho);
	for(i=0; i<Nfock; i++) {
		rho->data[i*Nfock + i].real = calc_thermal(i, nbr_photons);
		//printf("rho init %f\n",rho->data[i*Nfock + i].real);
	}

	//printf("trace init %Lf\n",calc_trace(rho));

	/*for (i=0; i<Nfock*Nfock; i++)
		printf("%f \n",rho->data[i].real);
	*/

	matrix *proj_sum2 = make_matrix(Nfock);

	char filename[100];
	snprintf(filename,sizeof(filename),"projectors_N%li_xlim%f_spacing%f_gridpts%li.bin",Nfock,xp_lim_pos, grid_spacing, nbr_grid_pts);



	// calculate projectors
	/*		printf("calculating projectors\n");
	calc_projs2(nbr_grid_pts, Nfock, grid_spacing, x, p, rho,  proj_sum2, filename);
	*/
				printf("calculating projector sum from saved projectors\n");
	calc_sum(Nfock, nbr_grid_pts, proj_sum2, filename);



	// print projector sum diagonal elements to file (all should be 1)
		f = fopen("diagonals.dat","w");
	for (i=0; i<Nfock; i++)
		fprintf(f,"%f\n",proj_sum2->data[i*Nfock + i].real);
	fclose(f);

	//	set_zero(rho);
	// read old rho to use as initial for maxlik
		/*	f = fopen("reconst_state_N500_maxlik2100_init1.dat", "r");
		for (i=0; i<Nfock*Nfock; i++){
			fscanf(f,"%lf %lf",&rho->data[i].real,&rho->data[i].imag);
			}
			fclose(f); */
	//print to check
	/*	for (i=0; i<10; i++){
			printf("rho %lf %lf\n",rho->data[i].real,rho->data[i].imag);
			}*/


	// set initial rho for maxlik
			set_zero(rho);
	for(i=0; i<Nfock; i++) {
		rho->data[i*Nfock + i].real = 1.0/Nfock;
		}
    printf("trace rho init %Lf\n",calc_trace(rho));

	maximum_likelihood(maxlik_iterations, cutoff, nbr_grid_pts, Nfock, H, rho, filename,proj_sum2);


	// print entire state
	f = fopen("reconst_state.dat", "w");
	fprintf(f,"Nfock %li\n",Nfock);
	for (i=0; i<Nfock*Nfock; i++){
		fprintf(f,"%f %f\n",rho->data[i].real,rho->data[i].imag);
	}
	fclose(f);

	// calculate trace (trace/Nfock should be 1)
	trace = calc_trace(proj_sum2);
   printf("trace state %Lf\n",trace/Nfock);

	// testprint sum
	//printf("1st element of sum %f\n",proj_sum2->data[0].real);

	// print projector sum diagonal elements to file (all should be 1)
	f = fopen("diagonals_state.dat","w");
		for (i=0; i<Nfock; i++)
		fprintf(f,"%f\n",proj_sum2->data[i*Nfock + i].real);
	fclose(f);


	free(H); free(x); free(p);

	mkl_free(rho->data); mkl_free(rho);
	//	mkl_free(proj_sum2->data); mkl_free(proj_sum2);
	/*
	for (i=0; i<nbr_of_projs; i++){
		mkl_free(Pi_diag[i]);
	}
	mkl_free(Pi_diag);*/
}
