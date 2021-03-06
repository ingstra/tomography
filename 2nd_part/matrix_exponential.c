# include <complex.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "matrix_exponential.h"
# include "c8lib.h"
# include "r8lib.h"
# include "functions.h"

/******************************************************************************/

void expm(matrix *A, matrix *e)

/******************************************************************************/
/*
  Parameters:

  Input, int N, the dimension of the matrix.

  Input, double complex A[N*N], the matrix.

  Output, e
*/
{
	int n = A->rows;


	matrix *A2 = make_matrix(n);
	matrix *X = make_matrix(n);

	double a_norm;
	double c;
	double e_tmp;
	double alpha, beta;
	alpha = 1; beta = 0;
	int ee;
	int k;
	const double one = 1.0;
	int p;
	const int q = 6;
	int s;
	int i;
	double t;


	//a2 = c8mat_copy_new ( n, n, a );
	cblas_zcopy(n*n, A->data, 1, A2->data, 1);

	//a_norm = c8mat_norm_li ( n, n, a2 );
	a_norm = LAPACKE_zlange(LAPACK_ROW_MAJOR, 'I', n, n, A2->data, n);
	//	printf("a_norm %f\n",a_norm);

	//ee = ( int ) ( r8_log_2 ( a_norm ) ) + 1;
	vdLog2(1,&a_norm,&e_tmp);
	ee = (int)e_tmp + 1;

	s = i4_max ( 0, ee + 1 );

	t = 1.0 / pow ( 2.0, s );

	//c8mat_scale_r8 ( n, n, t, a2 );
	cblas_zdscal(n*n, t, A2->data, 1);

	//printf("A2 0 %f\n",A2->data[0].real);

	//x = c8mat_copy_new ( n, n, a2 );
	cblas_zcopy(n*n, A2->data, 1, X->data, 1);

	//printf("X 1st 0 %f\n",A2->data[0].real);

	c = 0.5;
	//e = c8mat_identity_new ( n );
	for (i=0; i<n; i++)
		e->data[i*n + i].real = 1.0;

	//c8mat_add_r8 ( n, n, one, e, c, a2, e );
	matrix *tmp_a2=make_matrix(n);
	cblas_zcopy(n*n, A2->data, 1, tmp_a2->data, 1);
	matrix *tmp_e = make_matrix(n);
	cblas_zcopy(n*n, e->data, 1, tmp_e->data, 1);

	cblas_zdscal(n*n, one, tmp_e->data, 1);
	cblas_zdscal(n*n, c, tmp_a2->data, 1);
	vzAdd(n*n, tmp_a2->data, tmp_e->data, e->data);

	//printf("e 0 %f\n",e->data[0].real);

	//d = c8mat_identity_new ( n );
	matrix *d = make_matrix(n);
	for (i=0; i<n; i++)
		d->data[i*n + i].real = 1.0;


	//c8mat_add_r8 ( n, n, one, d, -c, a2, d );;
	cblas_zcopy(n*n, A2->data, 1, tmp_a2->data, 1);
	matrix *tmp_d = make_matrix(n);
	cblas_zcopy(n*n, d->data, 1, tmp_d->data, 1);

	cblas_zdscal(n*n, one, tmp_d->data, 1);
	cblas_zdscal(n*n, -c, tmp_a2->data, 1);
	vzAdd(n*n, tmp_a2->data, tmp_d->data, d->data);

	//printf("d 0 %f\n",d->data[0].real);

	matrix *tmp_x = make_matrix(n);
	p = 1;

	for ( k = 2; k <= q; k++ )
		{
			c = c * ( double ) ( q - k + 1 ) / ( double ) ( k * ( 2 * q - k + 1 ) );

			//c8mat_mm ( n, n, n, a2, x, x );
			set_zero(tmp_x);
			cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, A2->data, n, X->data, n, &beta, tmp_x->data, n);
			//	printf("X  %f\n",X->data[0].real);
			//printf("A2  %f\n",A2->data[0].real);
			//printf("tmp_x 1111 %f\n",tmp_x->data[0].real);
			cblas_zcopy(n*n, tmp_x->data, 1, X->data, 1);
			cblas_zcopy(n*n, e->data, 1, tmp_e->data, 1);
			//	printf("tmp_x 1st 0 %f\n",tmp_x->data[0].real);
			//c8mat_add_r8 ( n, n, c, x, one, e, e );
			cblas_zdscal(n*n, c, tmp_x->data, 1);
			cblas_zdscal(n*n, one, tmp_e->data, 1);
			vzAdd(n*n, tmp_x->data, tmp_e->data, e->data);
			/*	printf("tmp_x 0 %f\n",tmp_x->data[0].real);
				printf("tmp_e 0 %f\n",tmp_e->data[0].real);
				printf("e again 0 %f\n",e->data[0].real);
			*/
			if ( p )
				{
					//c8mat_add_r8 ( n, n, c, x, one, d, d );
					cblas_zcopy(n*n, X->data, 1, tmp_x->data, 1);
					cblas_zcopy(n*n, d->data, 1, tmp_d->data, 1);
					cblas_zdscal(n*n, c, tmp_x->data, 1);
					cblas_zdscal(n*n, one, tmp_d->data, 1);
					vzAdd(n*n, tmp_x->data, tmp_d->data, d->data);
					//printf("d 2 %f\n",d->data[0].real);

				}
			else
				{
					//c8mat_add_r8 ( n, n, -c, x, one, d, d );
					cblas_zcopy(n*n, X->data, 1, tmp_x->data, 1);
					cblas_zcopy(n*n, d->data, 1, tmp_d->data, 1);
					cblas_zdscal(n*n, -c, tmp_x->data, 1);
					cblas_zdscal(n*n, one, tmp_d->data, 1);
					vzAdd(n*n, tmp_x->data, tmp_d->data, d->data);
					//printf("d more %f\n",d->data[0].real);

				}

			p = !p;
		}
	/*
	  E -> inverse(D) * E
	*/
	//	printf("e 1st  %f\n",e->data[0].real);
	//c8mat_minvm ( n, n, d, e, e );
	int *ipiv = calloc(n*n,sizeof(int));
	matrix *tmp_lu = make_matrix(n);
	cblas_zcopy(n*n, d->data, 1, tmp_lu->data, 1);
	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, tmp_lu->data, n,ipiv);
	LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, tmp_lu->data, n, ipiv); // now tmp_lu is the inverse
	set_zero(tmp_e);
	cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, tmp_lu->data, n, e->data, n, &beta, tmp_e->data, n);
	cblas_zcopy(n*n, tmp_e->data, 1, e->data, 1);
	//	printf("e 2nd  %f\n",e->data[0].real);

	/*
	  E -> E^(2*S)
	*/
	for ( k = 1; k <= s; k++ )
		{
			//c8mat_mm ( n, n, n, e, e, e );

			/*	printf("k %i\n",k);
				printf("tmp_e before %f\n",tmp_e->data[0].real);
				printf("e before  %f\n",e->data[0].real);*/
			set_zero(e);
			cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, tmp_e->data, n, tmp_e->data, n, &beta, e->data, n);

			//			printf("tmp_e mid %f\n",tmp_e->data[0].real);
			//			printf("e mid  %f\n",e->data[0].real);
			cblas_zcopy(n*n, e->data, 1, tmp_e->data, 1);
			//	printf("tmp_e 0 %f\n",tmp_e->data[0].real);
			//printf("e more  %f\n",e->data[0].real);
		}

	// e is the result

		free(ipiv);

	mkl_free(A2->data);
	mkl_free(X->data);
	mkl_free(tmp_a2->data);
	mkl_free(tmp_e->data);
	mkl_free(tmp_d->data);
	mkl_free(d->data);
	mkl_free(tmp_lu->data);
	mkl_free(tmp_x->data);

	mkl_free(A2);
	mkl_free(X);
	/*mkl_free(tmp_a2);
	mkl_free(tmp_e);
	mkl_free(tmp_d);
	mkl_free(d);
	mkl_free(tmp_lu);
	mkl_free(tmp_x);

	*/
}

/******************************************************************************/

double complex *c8mat_expm1 ( int n, double complex a[] )

/******************************************************************************/
/*
  Purpose:

  C8MAT_EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.

  Licensing:

  This code is distributed under the GNU LGPL license.

  Modified:

  05 March 2013

  Author:

  Cleve Moler, Charles Van Loan

  Reference:

  Cleve Moler, Charles VanLoan,
  Nineteen Dubious Ways to Compute the Exponential of a Matrix,
  Twenty-Five Years Later,
  SIAM Review,
  Volume 45, Number 1, March 2003, pages 3-49.

  Parameters:

  Input, int N, the dimension of the matrix.

  Input, double complex A[N*N], the matrix.

  Output, double complex C8MAT_EXPM1[N*N], the estimate for exp(A).
*/
{
	double complex *a2;
	double a_norm;
	double c;
	double complex *d;
	double complex *e;
	int ee;
	int k;
	const double one = 1.0;
	int p;
	const int q = 6;
	int s;
	double t;
	double complex *x;

	a2 = c8mat_copy_new ( n, n, a );

	a_norm = c8mat_norm_li ( n, n, a2 );

	ee = ( int ) ( r8_log_2 ( a_norm ) ) + 1;

	s = i4_max ( 0, ee + 1 );

	t = 1.0 / pow ( 2.0, s );

	c8mat_scale_r8 ( n, n, t, a2 );
	//printf("a2 %f\n",creal(a2[0]));

	x = c8mat_copy_new ( n, n, a2 );
	//printf("x 1st ok %f\n",creal(x[0]));
	c = 0.5;

	e = c8mat_identity_new ( n );


	c8mat_add_r8 ( n, n, one, e, c, a2, e );

	// printf("e ok %f\n",creal(e[0]));

	d = c8mat_identity_new ( n );

	c8mat_add_r8 ( n, n, one, d, -c, a2, d );
	//printf("d ok %f\n",creal(d[0]));

	p = 1;

	for ( k = 2; k <= q; k++ )
		{
			c = c * ( double ) ( q - k + 1 ) / ( double ) ( k * ( 2 * q - k + 1 ) );

			c8mat_mm ( n, n, n, a2, x, x );
			//printf("x  %f\n",creal(x[0]));
			c8mat_add_r8 ( n, n, c, x, one, e, e );
			//printf("e2  %f\n",creal(e[0]));

			if ( p )
				{
					c8mat_add_r8 ( n, n, c, x, one, d, d );
					//printf("d  %f\n",creal(d[0]));
				}
			else
				{
					c8mat_add_r8 ( n, n, -c, x, one, d, d );
					//printf("d  %f\n",creal(d[0]));
				}

			p = !p;
		}
	/*
	  E -> inverse(D) * E
	*/
	c8mat_minvm ( n, n, d, e, e );
	// printf("e1 %f\n",creal(e[0]));

	/*
	  E -> E^(2*S)
	*/
	for ( k = 1; k <= s; k++ )
		{
			c8mat_mm ( n, n, n, e, e, e );
			//printf("e2 %f\n",creal(e[0]));
		}

	free ( a2 );
	free ( d );
	free ( x );

	return e;
}
