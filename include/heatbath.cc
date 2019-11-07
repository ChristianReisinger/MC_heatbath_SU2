#include <set>

#include <helper_functions.hh>

#include <geometry.hh>
#include <linear_algebra.hh>
#include <ranlux.hh>

#include <heatbath.hh>

void do_sweep(double* const gauge_field, int T, int L, double beta, const std::set<int>& fixed_timeslices) {
	for (int it = 0; it < T; it++) {
		bool on_boundary = fixed_timeslices.count(it);
		for (int ix = 0; ix < L; ix++) {
			for (int iy = 0; iy < L; iy++) {
				for (int iz = 0; iz < L; iz++) {
					for (int mu = 0; mu < 4; mu++) {

						if (on_boundary && mu != 0)
							continue;

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
						if (abs(norm - 1.0) > pow(10.0, -12.0))
							throw "Error: New link is not in SU(2).\n";
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
}
