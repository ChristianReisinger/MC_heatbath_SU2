// MC_heatbath.cc

// I do not remember, where I got these numbers --> be cautious!
// beta = 2.3   -->   a = 0.171 fm
// beta = 2.4   -->   a = 0.128 fm
// beta = 2.5   -->   a = 0.089 fm

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include <fields.hh>
#include <geometry.hh>
#include <io.hh>
#include <linear_algebra.hh>
#include <ranlux.hh>
#include <Wilson_loops.hh>

// lattice extensions
const int T, L;

double* gauge_field;

int main(int argc, char **argv) {
	using namespace std;

	if (argc != 9) {
		cout << "Usage: " << argv[0] << " <path> <beta> <T> <L> <seed> <cold/hot> <num_MC_sweeps_max> <num_MC_sweeps_out>\n";
		return 0;
	}

	string path(argv[1]);
	double beta;
	int seed, num_MC_sweeps_max, num_MC_sweeps_out;
	bool hot_start = true;

	stringstream args_ss;
	args_ss << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[7] << " " << argv[8];
	args_ss >> beta >> T >> L >> seed >> num_MC_sweeps_max >> num_MC_sweeps_out;

	if (string(argv[6]) == string("cold")) {
		hot_start = false;
	} else if (string(argv[6]) != string("hot")) {
		cerr << "Error: invalid parameter <cold/hot>\n";
		return 0;
	}

	if (beta <= 0.0) {
		cerr << "Error: beta <= 0.0\n";
		return 0;
	}

	if (T < 2 || L < 2) {
		cerr << "Error: T < 2 or L < 2\n";
		return 0;
	}

	if (seed < 1) {
		cerr << "Error: seed < 1\n";
		return 0;
	}

	if (num_MC_sweeps_max < 1 || num_MC_sweeps_out < 1) {
		cerr << "Error: num_MC_sweeps_max < 1 or num_MC_sweeps_out < 1.\n";
		return 0;
	}

	int config_id_digits = floor(log10(num_MC_sweeps_max)) + 1;

	InitializeRand(seed);

	Gauge_Field_Alloc(&gauge_field, T, L);

	if (hot_start) {
		cout << "Initializing hot start ...\n";
		Gauge_Field_Random(gauge_field, T, L);
	} else {
		cout << "Initializing cold start ...\n";
		Gauge_Field_Unity(gauge_field, T, L);
	}

	double P = Average_Plaquette(gauge_field, T, L);
	cout << "  ... <P> = " << P << ".\n";

	// Update links using the heatbath algorithm.

	ostringstream plaq_filename_oss;
	plaq_filename_oss << path << "/P.dat";
	ofstream plaquette_ofs(plaq_filename_oss.str());
	if (plaquette_ofs.fail()) {
		cerr << "Error: Unable to open average plaquette file.\n";
		return 0;
	}

	for (int MC_sweep_i = 1; MC_sweep_i <= num_MC_sweeps_max; MC_sweep_i++) {
		for (int it = 0; it < T; it++) {
			for (int ix = 0; ix < L; ix++) {
				for (int iy = 0; iy < L; iy++) {
					for (int iz = 0; iz < L; iz++) {
						for (int mu = 0; mu < 4; mu++) {

							if (ix == 0 && iy == 0 && iz == 0 && mu == 0) {
								cout << "sweep = " << setw(config_id_digits) << MC_sweep_i << ", it = " << setw(2) << it
										<< ", ix = " << setw(2) << ix << ", iy = " << setw(2) << iy
										<< ", iz = " << setw(2) << iz << ", mu = " << mu << " ...\n";
							}

// Cf. Montvay/Münster, sec. 7.3.1, page 399.

// S_l = "sum over staples", more precisely S = -\beta/N Re(Tr(U_l S_l)).

							double S_l[8];
							cm_eq_zero(S_l);

							double U1[8], U2[8];
							for (int nu = 0; nu < 4; nu++) {
								if (nu == mu)
									continue;

								int n[4] = { it, ix, iy, iz };
								int index1, index2, index3;

								// negative mu-nu staple
								n[nu] -= 1;
								index1 = ggi(get_index(n[0], n[1], n[2], n[3], T, L), nu);
								index2 = ggi(get_index(n[0], n[1], n[2], n[3], T, L), mu);
								n[mu] += 1;
								index3 = ggi(get_index(n[0], n[1], n[2], n[3], T, L), nu);
								cm_eq_cm_ti_cm(U1,
										gauge_field + index2,
										gauge_field + index3);
								cm_eq_cm_dag_ti_cm(U2, gauge_field + index1, U1);
								cm_pl_eq_cm(S_l, U2);

								// positive mu-nu staple
								n[0] = it;
								n[1] = ix;
								n[2] = iy;
								n[3] = iz;
								index1 = ggi(get_index(n[0], n[1], n[2], n[3], T, L), nu);
								n[nu] += 1;
								index2 = ggi(get_index(n[0], n[1], n[2], n[3], T, L), mu);
								n[mu] += 1;
								n[nu] -= 1;
								index3 = ggi(get_index(n[0], n[1], n[2], n[3], T, L), nu);
								cm_eq_cm_ti_cm_dag(U1,
										gauge_field + index2,
										gauge_field + index3);
								cm_eq_cm_ti_cm(U2, gauge_field + index1, U1);
								cm_pl_eq_cm(S_l, U2);
							}

							cm_dag_eq_cm(S_l);

// Cf. eq. (7.106).
// Step i) and ii) on page 401.
							double k = sqrt(S_l[0] * S_l[6] - S_l[1] * S_l[7] - S_l[2] * S_l[4] + S_l[3] * S_l[5]);
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

// Step iv) on page 401.
// U_0 = k S_l^{-1} = S_l^\dagger / k   (cf. eq. (7.105), (7.110))
							double U_0[8];
							cm_eq_cm_dag(U_0, S_l);
							cm_ti_eq_re(U_0, 1.0 / k);

// U_0l = a_0 + i \sigma_j a_j   (cf. eq. (7.110))
							double U_0l[8];
							cm_from_h(U_0l, a);
							cm_eq_cm_ti_cm(U1, U_0l, U_0);

// Project the link to SU(2) (roundoff errors accumulate during simulation).
							double h[4];
							h_from_cm(h, U1);
							norm = 1.0 / sqrt(pow(h[0], 2.0) + pow(h[1], 2.0) + pow(h[2], 2.0) + pow(h[3], 2.0));
							if (abs(norm - 1.0) > pow(10.0, -12.0)) {
								cerr << "Error: New link is not in SU(2).\n";
								return 1;
							}
							h[0] *= norm;
							h[1] *= norm;
							h[2] *= norm;
							h[3] *= norm;
							cm_from_h(U1, h);

							cm_eq_cm(gauge_field + ggi(get_index(it, ix, iy, iz, T, L), mu), U1);
						}
					}
				}
			}
		}

		cout << "  ... sweep " << MC_sweep_i << " - end ...\n";

		P = Average_Plaquette(gauge_field, T, L);
		cout << "  ... <P> = " << P << ".\n";
		plaquette_ofs << setw(5) << MC_sweep_i << " " << showpos << setprecision(5) << scientific << P << noshowpos << endl;

		if (MC_sweep_i % num_MC_sweeps_out == 0) {
			ostringstream config_filename_oss;
			config_filename_oss << path << "/conf." << setfill(0) << setw(config_id_digits) << MC_sweep_i;

			ostringstream header_oss(argv[0]);
			for (int arg_i = 1; arg_i < argc; ++arg_i)
				header_oss << " " << argv[arg_i];
			write_gauge_field(gauge_field, config_filename_oss.str().c_str(), T, L, header_oss.str().c_str());
		}
	}

	Gauge_Field_Free(&gauge_field);

	return EXIT_SUCCESS;
}

// ********************
