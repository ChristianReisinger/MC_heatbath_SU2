#include "smearing_techniques.hh"

#include "geometry.hh"
#include "geometry2.hh"
#include "linear_algebra.hh"

#include <global_defs.hh>

void APE_Smearing_Step(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha) {
	using namespace de_uni_frankfurt_itp::reisinger::latticetools_0719;

	double *unsmeared_gauge_field;
	Gauge_Field_Alloc(unsmeared_gauge_field, T, L);
	Gauge_Field_Copy(unsmeared_gauge_field, smeared_gauge_field, T, L);

	double U[SUN_elems], V[SUN_elems];

	for (int it = 0; it < T; ++it)
		for (int ix = 0; ix < L; ++ix)
			for (int iy = 0; iy < L; ++iy)
				for (int iz = 0; iz < L; ++iz)
					for (int mu = 1; mu < 4; ++mu) {
						double* U_new = smeared_gauge_field + ggi_n(it, ix, iy, iz, 4, mu, T, L);
						cm_eq_zero(U_new);
						for (int staple_dir = 1; staple_dir < 4; ++staple_dir) {
							if (staple_dir == mu)
								continue;

							//pos direction
							int n[4] = { it, ix, iy, iz };
							cm_eq_cm(U, unsmeared_gauge_field + ggi_n(n, 4, staple_dir, T, L));
							++n[staple_dir];
							cm_eq_cm_ti_cm(V, U, unsmeared_gauge_field + ggi_n(n, 4, mu, T, L));
							++n[mu];
							--n[staple_dir];
							cm_eq_cm_ti_cm_dag(U, V, unsmeared_gauge_field + ggi_n(n, 4, staple_dir, T, L));
							cm_pl_eq_cm(U_new, U);

							//neg direction
							--n[mu];
							--n[staple_dir];
							cm_eq_cm_dag(U, unsmeared_gauge_field + ggi_n(n, 4, staple_dir, T, L));
							cm_eq_cm_ti_cm(V, U, unsmeared_gauge_field + ggi_n(n, 4, mu, T, L));
							++n[mu];
							cm_eq_cm_ti_cm(U, V, unsmeared_gauge_field + ggi_n(n, 4, staple_dir, T, L));
							cm_pl_eq_cm(U_new, U);

							cm_ti_eq_re(U_new, APE_smearing_alpha);
						}
						cm_pl_eq_cm(U_new, unsmeared_gauge_field + ggi_n(it, ix, iy, iz, 4, mu, T, L));

						//Analytically, U_new is of the form (a, b; -b*, a*), but not numerically. Remove numerical errors:
						U_new[4] = -U_new[2];
						U_new[5] = U_new[3];
						U_new[6] = U_new[0];
						U_new[7] = -U_new[1];
						//Tested: for U_new  = (a, b; c, d), using (a+d*, b-c*; c-b*, d+a*)/2 gives the same results as using (a, b; -b*, a*)

						cm_proj(U_new);
					}

	Gauge_Field_Free(unsmeared_gauge_field);
}
