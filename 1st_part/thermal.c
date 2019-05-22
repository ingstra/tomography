#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>


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
	return result;
}


int main() {

	unsigned long int nbr_photons, m, N;

	nbr_photons=35;
	N=300;

	double *termal = malloc(N*sizeof(double));

	for (m=0; m<N; m++) {
		termal[m] = calc_thermal(m, nbr_photons);

	}

	FILE *f;
	f = fopen("thermal.dat","w");

	for (m=0; m<N; m++) {
		fprintf(f, "%f\n",termal[m]);

	}
	fclose(f);




	return 0;
}
