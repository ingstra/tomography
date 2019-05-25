#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "matrix.h"
#include <mpc.h>
#include <gmp.h>

long double complex assign_proj_element(double grid_spacing, double xi, double pj, unsigned long int m, unsigned long int n) {
    long double complex integration_approx = grid_spacing*grid_spacing;
    long double complex proj_element;

    /*** calculate factorials ***/
    // declare and initialize mpz integer factorial variables
    mpz_t factorial_m_z, factorial_n_z,factorial;
    mpz_init(factorial_m_z); mpz_init(factorial_n_z);mpz_init(factorial);
    // calculate integer factorials
    mpz_fac_ui(factorial_m_z,m);
    mpz_fac_ui(factorial_n_z,n);
	// multiply factorials
	mpz_mul(factorial,factorial_n_z,factorial_m_z);
	//initialize float, and convert integer factorial to float
	mpf_t factorial_f; mpf_init(factorial_f);
	mpf_set_z(factorial_f,factorial);
	// square root of factorials
	mpf_sqrt(factorial_f,factorial_f);


    // declare and initialize mpc complex factorial variables
    mpc_t factorial_m_c, factorial_n_c;
    mpc_init2(factorial_m_c,256); mpc_init2(factorial_n_c,256);
    // convert factorial from float to complex
    mpc_set_z(factorial_m_c, factorial_m_z, MPC_RNDNN);
    mpc_set_z(factorial_n_c, factorial_n_z, MPC_RNDNN);
	// take the square root of factorials
    mpc_sqrt(factorial_m_c,factorial_m_c,MPC_RNDNN);
    mpc_sqrt(factorial_n_c,factorial_n_c,MPC_RNDNN);

    /*** calculate powers ***/
    // declare and initialize mpc phase space point alpha from xi and pj
    mpc_t alpha, alpha_conj; mpc_init2(alpha,256); mpc_init2(alpha_conj,256);
    // assign value to alpha
    mpc_set_d_d(alpha,xi,pj,MPC_RNDNN);
    mpc_conj(alpha_conj,alpha,MPC_RNDNN);
    // declare and initialize variables to hold powers
    mpc_t P_m, P_n; mpc_init2(P_m,256); mpc_init2(P_n,256);
    // set powers
    mpc_pow_ui(P_m,alpha,m,MPC_RNDNN);
    mpc_pow_ui(P_n,alpha_conj,n,MPC_RNDNN);

	/*** exponential function part ***/
	// declare and initialize
	mpc_t exp,Xi_squared,Pj_squared,arg; mpc_init2(exp,256);mpc_init2(Xi_squared,256); mpc_init2(Pj_squared,256); mpc_init2(arg,256);
	// convert to mpc to be able to use mpc exponential function
	mpc_set_d(Xi_squared,-xi*xi,MPC_RNDNN);
	mpc_set_d(Pj_squared,-pj*pj,MPC_RNDNN);
	//calculate
	mpc_add(arg,Xi_squared,Pj_squared,MPC_RNDNN);
	mpc_exp(exp,arg,MPC_RNDNN);


    /*** put together some parts ***/
    // declare and initialize
    mpc_t part1, part2,total; mpc_init2(part1,256); mpc_init2(part2,256); mpc_init2(total,256);
    // put together
    mpc_div(part1,P_m,factorial_m_c,MPC_RNDNN);
    mpc_div(part2,P_n,factorial_n_c,MPC_RNDNN);
    mpc_mul(total,part1,part2,MPC_RNDNN);
	mpc_mul(total,total,exp,MPC_RNDNN);

    /*** convert to long double complex, and clear mpc and mpz variables ***/


	proj_element = mpc_get_ldc(total,MPC_RNDNN);
    mpc_clear(total); mpc_clear(part1); mpc_clear(part2); mpc_clear(P_m); mpc_clear(P_n);
    mpc_clear(alpha); mpc_clear(alpha_conj); mpc_clear(factorial_m_c); mpc_clear(factorial_n_c);
    mpz_clear(factorial_m_z); mpz_clear(factorial_n_z);
	mpc_clear(exp); mpc_clear(Xi_squared); mpc_clear(Pj_squared); mpc_clear(arg);

    proj_element = proj_element*integration_approx/M_PI;
    //printf("res %e %f %f %li %li  \n ",cimag(proj_element),xi,pj,m,n);
	return proj_element;
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
    double xp_lim, grid_spacing;
	long double trace;
	unsigned long int m, n;

	Nfock = 100;

	// projector Pi
	//double *Pi = malloc(Nfock * Nfock * sizeof(double));

	xp_lim = 10.0;
	nbr_grid_pts = 30;
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
	printf("elem1 sum %e \n", creal(proj_sum->data[0]));

}
