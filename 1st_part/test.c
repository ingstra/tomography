#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <mpc.h>
#include <gmp.h>



int main() {

	double grid_spacing = 0.8;
	double xi, pj;
	int m, n;

	xi = 15;
	pj = 15;

	m = 200;
	n = 200;

	mpc_t n1, n2, n3, d1, d2, full, test1, test2;
	mpc_init2(n1,256);
	mpc_init2(n2,256);
	mpc_init2(n3,256);
	mpc_init2(d1,256);
	mpc_init2(d2,256);
	mpc_init2(full,256);
	mpc_init2(test1,256);
	mpc_init2(test2,256);

	mpz_t factorial;
	mpz_init(factorial);
	unsigned long int N;
	N = 4;
	mpz_fac_ui(factorial,N);
	gmp_printf("factorial %Zd\n",factorial);



					/*n1 = cexp(-xi*xi -pj*pj);
	n2 = cpow(xi + I*pj,m);
	n3 = cpow(xi - I*pj,n);
	d1 = tgamma(m+1)*1e-200;
	d2 = csqrt(tgamma(m+1)*tgamma(n+1));
	full = cexp(-xi*xi -pj*pj) * cpow(xi + I*pj,m)*cpow(xi - I*pj,n) / (M_PI * csqrt(tgamma(m+1)*tgamma(n+1)));

	test1 = cpow(xi + I*pj,m) / csqrt(tgamma(m+1));
	test2 = cpow(xi - I*pj,n) / csqrt(tgamma(n+1));

	printf("n1 %e\n",creal(n1));
	printf("n2 %e\n",creal(n2));
	printf("n3 %e\n",creal(n3));
	printf("d1 %e\n",creal(d1));
	printf("d2 %e\n",creal(d2));
	/*
	printf("full %e\n",creal(full));
	printf("test1 %e\n",creal(test1));
	printf("test2 %e\n",creal(test2));
	printf("test1*test2 %e\n",creal(test1*test2*n1));
	*/

}
