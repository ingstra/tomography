#include "functions.h"


// create (allocate) COMPLEX zero-matrix
matrix* make_matrix(int rows, int cols) {
	matrix *m = mkl_malloc(sizeof(matrix),64);
	m->rows = rows;
	m->cols = cols;
	m->data = mkl_calloc(rows*cols, sizeof(MKL_Complex16),64);
	return m;
}

// linspace function. Output: pointer *xp
void linspace(int xp_lim, int nbr_grid_pts, double grid_spacing, double *xp) {
	int i;

	xp[0] = -xp_lim + grid_spacing/2;
	for (i=1; i<nbr_grid_pts; i++) {
		xp[i] = xp[i-1] + grid_spacing;
	}
}

// trace function: returns: trace
long double calc_trace(matrix *M) {
	int i;
	long double trace = 0;
	int dim = M->rows;

	for (i=0; i < dim; i++) {
		trace += M->data[dim * i + i].real;
		//printf("fcn trace %Le %e\n",trace, M->data[dim * i + i]);
	}
	return trace;
}

// set matrix to zero
void set_zero(matrix *M) {
	int i;
	for (i=0; i < M->rows * M->cols; i++) {
		M->data[i].real = 0;
		M->data[i].imag = 0;
	}
}


// calculate all projectors. Output: array of matrices Pi, and matrix proj_sum
void calc_all_projs(unsigned long int nbr_grid_pts, unsigned long int Nfock, double grid_spacing, double *x, double *p, matrix **Pi, matrix *proj_sum) {
unsigned long int i, j, m, n;

for (i=0; i<nbr_grid_pts; i++) {
		for (j=0; j< nbr_grid_pts; j++) {
			// allocate 0 matrix
			Pi[i*nbr_grid_pts + j] = make_matrix(Nfock,Nfock);
			// parallel for loop
            #pragma omp parallel for private(m,n)
			for (m=0; m<Nfock; m++) {
				for (n=0; n<Nfock; n++) {
					//printf("main x %e p %e i %i j %i m %i n %i nbrPi %i element %i \n",x[i],p[j],i,j,m,n,i*nbr_grid_pts + j,m*Nfock + n);
					// calc projector matrix elements
					Pi[i*nbr_grid_pts + j]->data[m*Nfock + n] = calc_proj_element(grid_spacing, x[i], p[j], m, n);
					//	printf("i %i, j %i, m %i, n %i, Pi %e\n",i,j,m,n, Pi[i*nbr_grid_pts + j]->data[m*Nfock + n].real);

					// sum all projectors
					proj_sum->data[Nfock*m + n].real += Pi[i*nbr_grid_pts + j]->data[m*Nfock + n].real;
				}
			}
		}
	}
// printf("PI %e\n",Pi[50*nbr_grid_pts+37]->data[Nfock].real);

}


// calculate projector matrix elements
MKL_Complex16 calc_proj_element(double grid_spacing, double xi, double pj, unsigned long int m, unsigned long int n) {
    long double complex integration_approx = grid_spacing*grid_spacing;
    long double complex proj_element;
	MKL_Complex16 element;

    /*** calculate factorials: factorial_c ***/
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
    mpc_t factorial_c;  mpc_init2(factorial_c,128);
    // convert factorial from float to complex
    mpc_set_f(factorial_c, factorial_f, MPC_RNDNN);

    /*** calculate powers: powers ***/
    // declare and initialize mpc phase space point alpha from xi and pj
    mpc_t alpha, alpha_conj; mpc_init2(alpha,128); mpc_init2(alpha_conj,128);
    // calc value to alpha
    mpc_set_d_d(alpha,xi,pj,MPC_RNDNN);
    mpc_conj(alpha_conj,alpha,MPC_RNDNN);
    // declare and initialize variables to hold powers
    mpc_t P_m, P_n, powers; mpc_init2(P_m,128); mpc_init2(P_n,128); mpc_init2(powers,128);
     // set powers
    mpc_pow_ui(P_m,alpha,m,MPC_RNDNN);
    mpc_pow_ui(P_n,alpha_conj,n,MPC_RNDNN);
	mpc_mul(powers,P_m,P_n,MPC_RNDNN);

	/*** exponential function part: exp ***/
	// declare and initialize
	mpc_t exp,Xi_squared,Pj_squared,arg; mpc_init2(exp,128);mpc_init2(Xi_squared,128); mpc_init2(Pj_squared,128); mpc_init2(arg,128);
	// convert to mpc to be able to use mpc exponential function
	mpc_set_d(Xi_squared,-xi*xi,MPC_RNDNN);
	mpc_set_d(Pj_squared,-pj*pj,MPC_RNDNN);
	//calculate
	mpc_add(arg,Xi_squared,Pj_squared,MPC_RNDNN);
	mpc_exp(exp,arg,MPC_RNDNN);


    /*** put together parts ***/
    // declare and initialize
    mpc_t total;   mpc_init2(total,128);
    // put together
	mpc_div(total,powers,factorial_c,MPC_RNDNN);
	mpc_mul(total,total,exp,MPC_RNDNN);


    /*** convert to long double complex, and clear mpc and mpz variables ***/
	proj_element = mpc_get_ldc(total,MPC_RNDNN);

    mpc_clear(total); mpc_clear(factorial_c);  mpc_clear(P_m); mpc_clear(P_n);
    mpc_clear(alpha); mpc_clear(alpha_conj); mpc_clear(powers); mpf_clear(factorial_f);
    mpz_clear(factorial_m_z); mpz_clear(factorial_n_z); mpz_clear(factorial);
	mpc_clear(exp); mpc_clear(Xi_squared); mpc_clear(Pj_squared); mpc_clear(arg);

    proj_element = proj_element*integration_approx/M_PI;
    //printf("res %e %f %f %li %li  \n ",cimag(proj_element),xi,pj,m,n);
	element.real = creal(proj_element);
	element.imag = cimag(proj_element);
	return element;
}


/*
  Input:
         iterations: number of iterations
		 cutoff: smallest number pr can be to be used
		 nbr_grid_pts
		 Nfock
         **Pi : projectors
         *rho: initial state
  Output:
         *rho: final state

 */
void maximum_likelihood(int iterations,long double cutoff, unsigned long int nbr_grid_pts, unsigned long int Nfock,  matrix **Pi, matrix *rho) {

	matrix *pr_matrix = make_matrix(Nfock,Nfock);
	// for frobenius diff
	matrix *rho_prev = make_matrix(Nfock,Nfock);
	matrix *rho_diff = make_matrix(Nfock,Nfock);
	matrix *rho_diff_dagger = make_matrix(Nfock,Nfock);
	matrix *rho_diff_prod = make_matrix(Nfock,Nfock);

	matrix *R = make_matrix(Nfock,Nfock);
	matrix *tmp1 =  make_matrix(Nfock,Nfock);

	double alpha = 1;
	double beta = 0;
	int i, j, k;

	long double rho_trace, pr, frobenius;
	MKL_Complex16 a; a.real=1; a.imag=0;


	FILE *f;

	f = fopen("frobenius.dat","w");
	printf("Starting maxlik\n");


	/********** MAXLIK *****************/
	for (k=0; k<iterations; k++) {
		// zero R
		set_zero(R);

		for (i=0; i<nbr_grid_pts; i++) {
			for (j=0; j<nbr_grid_pts; j++) {

				set_zero(pr_matrix);
				cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, rho->data, Nfock, Pi[i*nbr_grid_pts + j]->data, Nfock, &beta, pr_matrix->data, Nfock);

				// calculate probability pr
				pr = calc_trace(pr_matrix);

				if (pr > cutoff ) {
					// cblas_zdscal updates the scaled vector, in this case Pi
					cblas_zdscal(Nfock*Nfock,  H[i*nbr_grid_pts + j] / pr, Pi[i*nbr_grid_pts + j]->data, 1);

					// add to R
				 	vzAdd(Nfock*Nfock, R->data, Pi[i*nbr_grid_pts + j]->data, R->data);

				} // end cutoff statement
			}
		}

		// save old rho by copying rho to rho_prev
		cblas_zcopy(Nfock*Nfock, rho->data, 1, rho_prev->data, 1);

		/*	for (q=0;q<9;q++) {
			printf("rho %f, rho_prev %f \n",rho->data[q].real, rho_prev->data[q].real);
			} */


		// update rho
		set_zero(tmp1);
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, rho->data, Nfock, R->data, Nfock, &beta, tmp1->data, Nfock);
		set_zero(rho);
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, R->data, Nfock, tmp1->data, Nfock, &beta, rho->data, Nfock);
		//normalize rho
		rho_trace=calc_trace(rho);
		cblas_zdscal(Nfock*Nfock, 1.0/rho_trace, rho->data, 1);



		/**** calculate frobenius norm of matrix difference  every tenth step  ****/
		if (k%10 == 0) {

		// calculate difference, store in rho_diff, and copy
		vzSub(Nfock*Nfock, rho_prev->data, rho->data, rho_diff->data);

		// calculate hermitian conjugate
		mkl_zomatcopy('R', 'C', Nfock, Nfock, a, rho_diff->data, Nfock, rho_diff_dagger->data, Nfock);

		// matrix multiplication
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, rho_diff->data, Nfock, rho_diff_dagger->data, Nfock, &beta, rho_diff_prod->data, Nfock);

		frobenius = calc_trace(rho_diff_prod);
		frobenius = sqrtl(frobenius);
		printf(" frobenius %Le\n", frobenius);
		fprintf(f,"%i %Le\n",k, frobenius);

		} // end frobenius


	} // end maxlik iterations loop
	fclose(f);



}
