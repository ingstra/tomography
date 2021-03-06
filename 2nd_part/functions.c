#include "functions.h"
#include "matrix_exponential.h"

// create (allocate) COMPLEX zero-matrix
matrix* make_matrix(int rows) {
	matrix *m = mkl_malloc(sizeof(matrix),64);
	m->rows = rows;
	m->data = mkl_calloc(rows*rows, sizeof(MKL_Complex16),64);
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

	//printf("trace fcn %Lf\n",trace);
	return trace;
}

// set matrix to zero
void set_zero(matrix *M) {
	int i;
	for (i=0; i < M->rows * M->rows; i++) {
		M->data[i].real = 0;
		M->data[i].imag = 0;
	}
}

matrix *creation_op(int N) {
	int i;
	matrix *a = make_matrix(N);
	for (i=0; i<N-1; i++) {
		a->data[i*N + i + 1].real = sqrt((double)(i+1));
	}
	return a;
}

/*
  returns (diagonal) matrix elements m of a thermal state
*/
double calc_thermal(unsigned long int m, unsigned long int nbr_photons) {
	double result;

	mpz_t numerator_z, denominator_z;
	mpz_init(numerator_z); mpz_init(denominator_z);
	mpf_t numerator_f, denominator_f, result_f;
	mpf_init(numerator_f); mpf_init(denominator_f);
	mpf_init(result_f);

	mpz_ui_pow_ui(numerator_z, nbr_photons, m);
	mpz_ui_pow_ui(denominator_z, 1 + nbr_photons, m+1);

	mpf_set_z(denominator_f, denominator_z);
	mpf_set_z(numerator_f, numerator_z);

	mpf_div(result_f, numerator_f, denominator_f);

	result = mpf_get_d(result_f);

	mpz_clear(numerator_z); mpz_clear(denominator_z);
	mpf_clear(numerator_f); mpf_clear(denominator_f); mpf_clear(result_f);

	return result;
}


/*
  D = expm(alpha*a - conj(alpha)*a_dagger)
*/
/*matrix *displacement(double x, double p, int N){
// creation/annihilation operators
matrix *D = make_matrix(N,N);
matrix *sub = make_matrix(N,N);
matrix *a = creation_op(N);
matrix *a_dagger = make_matrix(N,N);
MKL_Complex16 scale; scale.real=1; scale.imag=0; // hermitian conjugate, no scaling
mkl_zomatcopy('R', 'C', N, N, scale, a->data, N, a_dagger->data, N);

// alpha
MKL_Complex16 alpha; alpha.real = x; alpha.imag = p;
MKL_Complex16 alpha_conj; alpha_conj.real = x; alpha_conj.imag = -p;

// multiply (scale)
mkl_zimatcopy('R', 'N', N, N, alpha_conj, a->data, N, N);
mkl_zimatcopy('R', 'N', N, N, alpha, a_dagger->data, N, N);
//subtract
vzSub(N*N, a_dagger->data, a->data, sub->data);

D = expm(sub);
mkl_free(a->data); mkl_free(a);
mkl_free(a_dagger->data); mkl_free(a_dagger);
mkl_free(sub->data); mkl_free(sub);
return D;
}
*/

// test if the projectors sum to the identity, without storing all of them
void test_projs2(unsigned long int nbr_grid_pts, unsigned long int N, double grid_spacing, double *x, double *p, matrix *rho_thermal, matrix *proj_sum) {

	unsigned long int i, j, k;
	matrix *D = make_matrix(N);
	matrix *D_dagger = make_matrix(N);
	matrix *tmp1 =  make_matrix(N); // for intermediate step in multipication
	MKL_Complex16 scale; scale.real=1; scale.imag=0; // for hermitian conjugate, no scaling
	matrix *Pi_tmp =  make_matrix(N);
	double const1 = 1;
	double const2 = 0;
	set_zero(proj_sum);
	//printf("D1 %f\n",D->data[0]);
	for (i=0; i<nbr_grid_pts; i++) {
		for (j=0; j< nbr_grid_pts; j++) {

			//	printf("hej\n");
		    set_zero(D);
			displacement(x[i], p[j], N,D);
		   	//printf("xi %f pj %f i %i j %i \n",x[i],p[j],i, j);
				/*	printf("%f %f %f \n",D->data[k*3].real,D->data[k*3+1].real,D->data[k*3+2].real);
					printf("D above \n");
				*/
				//	printf("hej2 %f\n",D->data[0]);
				//	printf("hej2\n");
				mkl_zomatcopy('R', 'C', N, N, scale, D->data, N, D_dagger->data, N);
			//	printf("hej3\n");
			// multiplication
			set_zero(tmp1);
			cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &const1, rho_thermal->data, N, D_dagger->data, N, &const2, tmp1->data, N);

			/*   	for (k=0;k<N;k++)
					printf("%f %f %f \n",tmp1->data[k*3].real,tmp1->data[k*3+1].real,tmp1->data[k*3+2].real);
					printf("tmp1 above \n");
			*/
			set_zero(Pi_tmp);
			cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &const1, D->data, N, tmp1->data, N, &const2, Pi_tmp->data, N);

			/*	for (k=0;k<N;k++){
				printf("%f %f %f \n",Pi[i*nbr_grid_pts + j]->data[k*3].real,Pi[i*nbr_grid_pts + j]->data[k*3+1].real,Pi[i*nbr_grid_pts + j]->data[k*3+2].real);}
				printf("Pi above \n");
			*/
			// divide by pi, and multipy with grid_spacing^ (because of integral approx)
			cblas_zdscal(N*N, grid_spacing*grid_spacing/M_PI, Pi_tmp->data, 1);

			//printf("xi %f pj %f i %i j %i \n",x[i],p[j],i, j);
				/*	for (k=0;k<N;k++){
				printf("%f %f %f \n",Pi[i*nbr_grid_pts + j]->data[k*3].real,Pi[i*nbr_grid_pts + j]->data[k*3+1].real,Pi[i*nbr_grid_pts + j]->data[k*3+2].real);}
				printf("\n");
			*/
			//printf("%e\n",Pi[i*nbr_grid_pts + j]->data[0].real);

			// sum all projectors
			vzAdd(N*N, proj_sum->data, Pi_tmp->data, proj_sum->data);

		}
	}
	mkl_free(D->data); mkl_free(D);
	mkl_free(D_dagger->data); mkl_free(D_dagger);
	mkl_free(tmp1->data); mkl_free(tmp1);
	// printf("PI %e\n",Pi[50*nbr_grid_pts+37]->data[Nfock].real);

}
/* calculate all projectors for state reconstruction
   Output: array of matrices Pi
   matrix proj_sum
*/
void calc_projs2(unsigned long int nbr_grid_pts, unsigned long int N, double grid_spacing, double *x, double *p, matrix *rho_thermal, matrix *proj_sum, char *filename) {

	unsigned long int i, j, k;
	matrix *D = make_matrix(N);
	matrix *D_dagger = make_matrix(N);
	matrix *tmp1 =  make_matrix(N); // for intermediate step in multipication
	MKL_Complex16 scale; scale.real=1; scale.imag=0; // for hermitian conjugate, no scaling
	double const1 = 1;
	double const2 = 0;
	set_zero(proj_sum);
	FILE *f = fopen(filename,"wb");

	matrix *Pi_tmp = make_matrix(N);

	//printf("D1 %f\n",D->data[0]);
	for (i=0; i<nbr_grid_pts; i++) {
		for (j=0; j< nbr_grid_pts; j++) {

			// allocate 0 matrix
			//Pi[i*nbr_grid_pts + j] = make_matrix(N); old
			//	printf("hej\n");

		    set_zero(D);
			set_zero(Pi_tmp);
			displacement(x[i], p[j], N,D);

		   	//printf("xi %f pj %f i %i j %i \n",x[i],p[j],i, j);
		   	//for (k=0;k<N;k++)
				/*	printf("%f %f %f \n",D->data[k*3].real,D->data[k*3+1].real,D->data[k*3+2].real);
					printf("D above \n");
				*/
				//	printf("hej2 %f\n",D->data[0]);
				//	printf("hej2\n");

				mkl_zomatcopy('R', 'C', N, N, scale, D->data, N, D_dagger->data, N);

			//	printf("hej3\n");
			// multiplication

			set_zero(tmp1);
			cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &const1, rho_thermal->data, N, D_dagger->data, N, &const2, tmp1->data, N);

			/*   	for (k=0;k<N;k++)
					printf("%f %f %f \n",tmp1->data[k*3].real,tmp1->data[k*3+1].real,tmp1->data[k*3+2].real);
					printf("tmp1 above \n");
			*/

			cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &const1, D->data, N, tmp1->data, N, &const2, Pi_tmp->data, N);

			/*	for (k=0;k<N;k++){
				printf("%f %f %f \n",Pi[i*nbr_grid_pts + j]->data[k*3].real,Pi[i*nbr_grid_pts + j]->data[k*3+1].real,Pi[i*nbr_grid_pts + j]->data[k*3+2].real);}
				printf("Pi above \n");
			*/
			// divide by pi, and multipy with grid_spacing^ (because of integral approx)

			cblas_zdscal(N*N, grid_spacing*grid_spacing/M_PI, Pi_tmp->data, 1);

				// print to file
		    fwrite(Pi_tmp->data, sizeof(MKL_Complex16),N*N,f);
			printf("trace in calc %lf\n",calc_trace(Pi_tmp));

			//printf("xi %f pj %f i %i j %i \n",x[i],p[j],i, j);
				/*	for (k=0;k<N;k++){
				printf("%f %f %f \n",Pi[i*nbr_grid_pts + j]->data[k*3].real,Pi[i*nbr_grid_pts + j]->data[k*3+1].real,Pi[i*nbr_grid_pts + j]->data[k*3+2].real);}
				printf("\n");
			*/
			//printf("%e\n",Pi[i*nbr_grid_pts + j]->data[0].real);

			// sum all projectors
			vzAdd(N*N, proj_sum->data, Pi_tmp->data, proj_sum->data);

		}
	}
	mkl_free(D->data); mkl_free(D);
	mkl_free(D_dagger->data); mkl_free(D_dagger);
	mkl_free(tmp1->data); mkl_free(tmp1);
	// printf("PI %e\n",Pi[50*nbr_grid_pts+37]->data[Nfock].real);
	fclose(f);
}

void calc_sum(int Nfock, int nbr_grid_pts, matrix *proj_sum, char *filename) {

	unsigned long int i, j, info;
	FILE *proj_file;
	proj_file = fopen(filename,"rb");
	if (proj_file == NULL) {
		printf("open for read file error, filename %s\n",filename);
		exit(0);
	}

	matrix *Pi = make_matrix(Nfock);


	for (i=0; i<nbr_grid_pts; i++) {
		for (j=0; j< nbr_grid_pts; j++) {

			info = fread(Pi->data, sizeof(MKL_Complex16), Nfock*Nfock, proj_file);
			// sum all projectors
			vzAdd(Nfock*Nfock, proj_sum->data, Pi->data, proj_sum->data);
			//printf("%e\n",proj_sum->data[1]);

		}
	}

	fclose(proj_file);
}

/*
  D = expm(alpha*a - conj(alpha)*a_dagger)
*/
void displacement(double x, double p, int N, matrix *D){
	// creation/annihilation operators
	matrix *sub = make_matrix(N);
	matrix *a = creation_op(N);
	matrix *a_dagger = make_matrix(N);
	MKL_Complex16 scale; scale.real=1; scale.imag=0; // hermitian conjugate, no scaling
	mkl_zomatcopy('R', 'C', N, N, scale, a->data, N, a_dagger->data, N);

	// alpha
	MKL_Complex16 alpha; alpha.real = x; alpha.imag = p;
	MKL_Complex16 alpha_conj; alpha_conj.real = x; alpha_conj.imag = -p;

	// multiply (scale)
    mkl_zimatcopy('R', 'N', N, N, alpha_conj, a->data, N, N);
	mkl_zimatcopy('R', 'N', N, N, alpha, a_dagger->data, N, N);
	//subtract
	vzSub(N*N, a_dagger->data, a->data, sub->data);

	expm(sub,D);

	mkl_free(a->data); mkl_free(a);
	mkl_free(a_dagger->data); mkl_free(a_dagger);
	mkl_free(sub->data); mkl_free(sub);
}



// calculate projector matrix elements for noise reconstruction
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


/* calculate all projectors for noise reconstruction
   Output: array of matrices Pi
   matrix proj_sum
*/

void calc_all_projs(unsigned long int nbr_grid_pts, unsigned long int Nfock, double grid_spacing, double *x, double *p, matrix **Pi, matrix *proj_sum) {
	unsigned long int i, j, m, n;
	FILE *f;
	size_t info;
	f = fopen("projectors.bin","wb");

	for (i=0; i<nbr_grid_pts; i++) {
		for (j=0; j< nbr_grid_pts; j++) {
			// allocate 0 matrix
			Pi[i*nbr_grid_pts + j] = make_matrix(Nfock);
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

			// print to file
			info = fwrite(Pi[i*nbr_grid_pts + j]->data, sizeof(MKL_Complex16),Nfock*Nfock,f);
			//printf("info write %i\n",info);
		}
	}

	// printf("PI %e\n",Pi[50*nbr_grid_pts+37]->data[Nfock].real);
	fclose(f);
}

/*
  Input:
  iterations: number of iterations
  cutoff: smallest number pr can be to be used
  nbr_grid_pts
  Nfock
  **Pi : projectors
  *histogram: histogram data
  *rho: initial state
  Output:
  *rho: final state

  */
void maximum_likelihood(int iterations,long double cutoff, unsigned long int nbr_grid_pts, unsigned long int Nfock, double *histogram, matrix *rho, char *filename, matrix *proj_sum) {

	matrix *pr_matrix = make_matrix(Nfock);
	// for frobenius diff
	matrix *rho_prev = make_matrix(Nfock);
	matrix *rho_diff = make_matrix(Nfock);
	matrix *rho_diff_dagger = make_matrix(Nfock);
	matrix *rho_diff_prod = make_matrix(Nfock);

	matrix *R = make_matrix(Nfock);
	matrix *tmp1 =  make_matrix(Nfock);

	double alpha = 1;
	double beta = 0;
	int i, j, k, q;

	long double rho_trace, pr, frobenius;
	MKL_Complex16 a; a.real=1; a.imag=0;
	size_t info;

	FILE *f;

	f = fopen("frobenius.dat","w");
	printf("Starting maxlik\n");

	matrix *Pi_test = make_matrix(Nfock);

	FILE *proj_file;
	proj_file = fopen(filename,"rb");
	if (proj_file == NULL) {
		printf("open for read file error maxlik, filename %s\n",filename);
		exit(0);
	}

	/********** MAXLIK *****************/
	for (k=0; k<iterations; k++) {
		// zero R
		set_zero(R);
		rewind(proj_file);

		//printf("trace0 %Lf rho0 %f\n",calc_trace(rho),rho->data[0].real);

		for (i=0; i<nbr_grid_pts; i++) {

			for (j=0; j<nbr_grid_pts; j++) {

				// Read from file

				info = fread(Pi_test->data, sizeof(MKL_Complex16), Nfock*Nfock, proj_file);
				//printf("Number of projector elements read from file: %i\n",info);
				//printf("trace000 %Lf\n",calc_trace(rho));
				set_zero(pr_matrix);
				cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, rho->data, Nfock, Pi_test->data, Nfock, &beta, pr_matrix->data, Nfock);


				// calculate probability pr
				pr = calc_trace(pr_matrix);


				if (pr > cutoff ) {
					//printf("pr %Le\n",pr);
					// cblas_zdscal updates the scaled vector, in this case Pi
					cblas_zdscal(Nfock*Nfock,  histogram[i*nbr_grid_pts + j] / pr, Pi_test->data, 1);
					//printf("trace R1 %Lf\n",calc_trace(R));
					//printf("trace Pi %Lf\n",calc_trace(Pi_test));
					// add to R
				 	vzAdd(Nfock*Nfock, R->data, Pi_test->data, R->data);

					//printf("trace R2 %Lf\n",calc_trace(R));

				} // end cutoff statement
			}

		}

		// save old rho by copying rho to rho_prev
		cblas_zcopy(Nfock*Nfock, rho->data, 1, rho_prev->data, 1);

		//	printf("trace1 %Lf\n",calc_trace(rho));

		// update rho
		set_zero(tmp1);
		cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, rho->data, Nfock, R->data, Nfock, &beta, tmp1->data, Nfock);

		//printf("tracetmp %Lf\n",calc_trace(tmp1));
		//printf("trace.15 %Lf\n",calc_trace(rho));


		set_zero(rho);
		cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, R->data, Nfock, tmp1->data, Nfock, &beta, rho->data, Nfock);

		// correction if projectors don't sum to unity
		// first invert sum
		int *ipiv = calloc(Nfock*Nfock,sizeof(int));
		matrix *sum_inv = make_matrix(Nfock);
		cblas_zcopy(Nfock*Nfock, proj_sum->data, 1, sum_inv->data, 1);
		LAPACKE_zgetrf(LAPACK_ROW_MAJOR, Nfock, Nfock, sum_inv->data, Nfock,ipiv);
		LAPACKE_zgetri(LAPACK_ROW_MAJOR, Nfock, sum_inv->data, Nfock, ipiv); // now sum_inv is the inverse of proj_sum
		// update rho with correction
			set_zero(tmp1);
			cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, rho->data, Nfock, sum_inv->data, Nfock, &beta, tmp1->data, Nfock);

			//	printf("trace %Lf\n",calc_trace(rho));

			set_zero(rho);
		cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, sum_inv->data, Nfock, tmp1->data, Nfock, &beta, rho->data, Nfock);



		//normalize rho
		rho_trace=calc_trace(rho);
		printf("trace before norm %Lf\n",rho_trace);
		cblas_zdscal(Nfock*Nfock, 1.0/rho_trace, rho->data, 1);
		printf("trace after norm %Lf\n",calc_trace(rho));



		/**** calculate frobenius norm of matrix difference  every tenth step  ****/
		if (k%10 == 0) {
			printf("Iteration %i\n",k);
			// calculate difference, store in rho_diff, and copy
			vzSub(Nfock*Nfock, rho_prev->data, rho->data, rho_diff->data);
			printf("trace frob %Lf\n",calc_trace(rho));
			// calculate hermitian conjugate
			mkl_zomatcopy('R', 'C', Nfock, Nfock, a, rho_diff->data, Nfock, rho_diff_dagger->data, Nfock);

			// matrix multiplication
			set_zero(rho_diff_prod);
			cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nfock, Nfock, Nfock, &alpha, rho_diff->data, Nfock, rho_diff_dagger->data, Nfock, &beta, rho_diff_prod->data, Nfock);

			frobenius = calc_trace(rho_diff_prod);
			frobenius = sqrtl(frobenius);
			printf(" frobenius %Le\n", frobenius);
			fprintf(f,"%i %Le\n",k, frobenius);

		} // end frobenius


	} // end maxlik iterations loop
	fclose(f);
	fclose(proj_file);


	mkl_free(pr_matrix->data);
	mkl_free(rho_prev->data);
	mkl_free(rho_diff->data);
	mkl_free(rho_diff_dagger->data);
	mkl_free(rho_diff_prod->data);
	mkl_free(R->data);
	mkl_free(tmp1->data);

	/********** END MAXLIK *****************/

}
