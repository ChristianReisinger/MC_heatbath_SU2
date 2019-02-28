// MC_heatbath.cc

// I do not remember, where I got these numbers --> be cautious!
// beta = 2.3   -->   a = 0.171 fm
// beta = 2.4   -->   a = 0.128 fm
// beta = 2.5   -->   a = 0.089 fm

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fields.hh>
#include <geometry.hh>
#include <io.hh>
#include <linear_algebra.hh>
#include <ranlux.hh>
#include <Wilson_loops.hh>

// Temporal extension of the lattice.
int T;

// Spatial extension of the lattice.
int L;

// Gauge field.
double *gauge_field;

// ********************

int main(int argc, char **argv) {
	FILE *file1;
	int i1, i2;
	char string1[1000], string2[1000];
	double SU2_1[8], SU2_2[8];

	int it, ix, iy, iz;

	// **********

	if (argc != 9) {
		fprintf(stderr, "%s <path> <beta> <T> <L> <seed> <cold/hot> <num_MC_sweeps_max> <num_MC_sweeps_out>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	// ***

	char path[1000];
	strcpy(path, argv[1]);

	// ***

	double beta = atof(argv[2]);

	if (beta <= 0.0) {
		fprintf(stderr, "Error: int main(...\n");
		fprintf(stderr, "beta <= 0.0.\n");
		exit(EXIT_FAILURE);
	}

	// ***

	T = atoi(argv[3]);
	L = atoi(argv[4]);

	if (T < 2 || L < 2) {
		fprintf(stderr, "Error: int main(...\n");
		fprintf(stderr, "T < 2 || L < 2.\n");
		exit(EXIT_FAILURE);
	}

	// ***

	int seed = atoi(argv[5]);

	if (seed < 1) {
		fprintf(stderr, "Error: int main(...\n");
		fprintf(stderr, "seed < 1.\n");
		exit(EXIT_FAILURE);
	}

	InitializeRand(seed);

	// ***

	bool hot_start;

	if (strcmp(argv[6], "cold") == 0) {
		hot_start = false;
	} else if (strcmp(argv[6], "hot") == 0) {
		hot_start = true;
	} else {
		fprintf(stderr, "Error: int main(...\n");
		fprintf(stderr, "<cold/hot> != cold && <cold/hot> != hot.\n");
		exit(EXIT_FAILURE);
	}

	// ***

	int num_MC_sweeps_max = atoi(argv[7]);
	int num_MC_sweeps_out = atoi(argv[8]);

	if (num_MC_sweeps_max < 1 || num_MC_sweeps_out < 1) {
		fprintf(stderr, "Error: int main(...\n");
		fprintf(stderr, "num_MC_sweeps_max < 1 || num_MC_sweeps_out < 1.\n");
		exit(EXIT_FAILURE);
	}

	// **********

	// Allocate memory for the gauge link configuration.
	Gauge_Field_Alloc(&gauge_field, T, L);

	if (hot_start == false)
	// Cold start.
			{
		fprintf(stderr, "Initializing cold start ...\n");
		Gauge_Field_Unity(gauge_field, T, L);
	} else
	// Hot start.
	{
		fprintf(stderr, "Initializing hot start ...\n");
		Gauge_Field_Random(gauge_field, T, L);
	}

	double P = Average_Plaquette(gauge_field, T, L);
	fprintf(stderr, "  ... <P> = %+.6lf.\n", P);
	fprintf(stdout, "%5d %+.5e\n", 0, P);

	// **********

	// Update links using the heatbath algorithm.

	sprintf(string1, "%sP.dat", path);

	if ((file1 = fopen(string1, "w")) == NULL) {
		fprintf(stderr, "Error: int main(...\n");
		fprintf(stderr, "Unable to open average plaquette file.\n");
		exit(EXIT_FAILURE);
	}

	int num_MC_sweeps;

	for (num_MC_sweeps = 0; num_MC_sweeps < num_MC_sweeps_max; num_MC_sweeps++) {
		for (it = 0; it < T; it++) {
			for (ix = 0; ix < L; ix++) {
				for (iy = 0; iy < L; iy++) {
					for (iz = 0; iz < L; iz++) {
						for (i1 = 0; i1 < 4; i1++) {
							if (ix == 0 && iy == 0 && iz == 0 && i1 == 0) {
								fprintf(stderr, "sweep = %5d, it = %2d, ix = %2d, iy = %2d, iz = %2d, i1 = %d ...\n",
										num_MC_sweeps + 1, it, ix, iy, iz, i1);
							}

// Cf. Montvay/Münster, sec. 7.3.1, page 399.

// *****

// S_l = "sum over staples", more precisely S = -\beta/N Re(Tr(U_l S_l)).

							double S_l[8];
							cm_eq_zero(S_l);

							for (i2 = 0; i2 < 4; i2++) {
								if (i2 == i1)
									continue;

								int i[4];
								int index1, index2, index3;

								// negative i1-i2 staple
								i[0] = it;
								i[1] = ix;
								i[2] = iy;
								i[3] = iz;
								i[i2] -= 1;
								index1 = ggi(get_index(i[0], i[1], i[2], i[3], T, L), i2);
								index2 = ggi(get_index(i[0], i[1], i[2], i[3], T, L), i1);
								i[i1] += 1;
								index3 = ggi(get_index(i[0], i[1], i[2], i[3], T, L), i2);
								cm_eq_cm_ti_cm(SU2_1, gauge_field + index2, gauge_field + index3);
								cm_eq_cm_dag_ti_cm(SU2_2, gauge_field + index1, SU2_1);
								cm_pl_eq_cm(S_l, SU2_2);

								// positive i1-i2 staple
								i[0] = it;
								i[1] = ix;
								i[2] = iy;
								i[3] = iz;
								index1 = ggi(get_index(i[0], i[1], i[2], i[3], T, L), i2);
								i[i2] += 1;
								index2 = ggi(get_index(i[0], i[1], i[2], i[3], T, L), i1);
								i[i1] += 1;
								i[i2] -= 1;
								index3 = ggi(get_index(i[0], i[1], i[2], i[3], T, L), i2);
								cm_eq_cm_ti_cm_dag(SU2_1, gauge_field + index2, gauge_field + index3);
								cm_eq_cm_ti_cm(SU2_2, gauge_field + index1, SU2_1);
								cm_pl_eq_cm(S_l, SU2_2);
							}

							cm_dag_eq_cm(S_l);

// *****

// Cf. eq. (7.106).
							double k = sqrt(S_l[0] * S_l[6] - S_l[1] * S_l[7] - S_l[2] * S_l[4] + S_l[3] * S_l[5]);

// *****

// Step i) and ii) on page 401.

							double beta_k = beta * k;

							double y_min = exp(-beta_k);
							double y_max = exp(+beta_k);
							double y;

							double a[4];

							while (1) {
								y = y_min + (y_max - y_min) * DRand();

								a[0] = log(y) / beta_k;

								if (DRand() <= sqrt(1.0 - pow(a[0], 2.0)))
									break;
							}

// *****

// Step iii) on page 401.

							double norm;

							while (1) {
								a[1] = 2.0 * DRand() - 1.0;
								a[2] = 2.0 * DRand() - 1.0;
								a[3] = 2.0 * DRand() - 1.0;

								norm = pow(a[1], 2.0) + pow(a[2], 2.0) + pow(a[3], 2.0);

								if (norm >= 0.0000000001 && norm <= 1.0)
									break;
							}

							norm = sqrt((1.0 - pow(a[0], 2.0)) / norm);

							a[1] *= norm;
							a[2] *= norm;
							a[3] *= norm;

// *****

// Step iv) on page 401.

// U_0 = k S_l^{-1} = S_l^\dagger / k   (cf. eq. (7.105), (7.110))
							double U_0[8];
							cm_eq_cm_dag(U_0, S_l);
							cm_ti_eq_re(U_0, 1.0 / k);

// U_0l = a_0 + i \sigma_j a_j   (cf. eq. (7.110))
							double U_0l[8];
							cm_from_h(U_0l, a);

							cm_eq_cm_ti_cm(SU2_1, U_0l, U_0);

// Project the link to SU(2) (roundoff errors accumulate during simulation).
							double h[4];
							h_from_cm(h, SU2_1);
							norm = 1.0 / sqrt(pow(h[0], 2.0) + pow(h[1], 2.0) + pow(h[2], 2.0) + pow(h[3], 2.0));
							if (fabs(norm - 1.0) > pow(10.0, -12.0)) {
								fprintf(stderr, "Error: int main(...\n");
								fprintf(stderr, "New link is not in SU(2).\n");
								exit(EXIT_FAILURE);
							}
							h[0] *= norm;
							h[1] *= norm;
							h[2] *= norm;
							h[3] *= norm;
							cm_from_h(SU2_1, h);

							cm_eq_cm(gauge_field + ggi(get_index(it, ix, iy, iz, T, L), i1), SU2_1);
						}
					}
				}
			}
		}

		fprintf(stderr, "  ... sweep %5d - end ...\n", num_MC_sweeps + 1);

		P = Average_Plaquette(gauge_field, T, L);
		fprintf(stderr, "  ... <P> = %+.6lf.\n", P);
		fprintf(file1, "%5d %+.5e\n", num_MC_sweeps + 1, P);
		fflush(file1);

		if ((num_MC_sweeps + 1) % num_MC_sweeps_out == 0) {
			sprintf(string1, "%sconf.%04d", path, num_MC_sweeps + 1);
			sprintf(string2, "%s %s %s %s %s %s %s %s %s", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7],
					argv[8]);
			write_gauge_field(gauge_field, string1, T, L, string2);
		}
	}

	fclose(file1);

	// **********

	Gauge_Field_Free(&gauge_field);

	// **********

	return EXIT_SUCCESS;
}

// ********************
