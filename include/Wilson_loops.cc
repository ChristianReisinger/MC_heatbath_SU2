// Wilson_loops.cc

#include "./Wilson_loops.hh"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./fields.hh"
#include "./geometry.hh"

// ********************

// Computes a Wilson loop.

// (it,ix,iy,iz): starting point.
// (dir1,dir2): the two directions.
// (ext1, ext2): the two extensions.

complex Wilson_Loop(double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dir1, int dir2, int ext1, int ext2) {
	complex c1;
	int i1;

	if (dir1 < 0 || dir1 > 3 || dir2 < 0 || dir2 > 3 || dir1 == dir2) {
		fprintf(stderr, "Error: complex Wilson_Loop(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8];

	cm_eq_id(U);

	// +dir_1-direction.

	for (i1 = 0; i1 < ext1; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), dir1));

		cm_eq_cm(U, U_tmp);

		if (dir1 == 0)
			it++;
		if (dir1 == 1)
			ix++;
		if (dir1 == 2)
			iy++;
		if (dir1 == 3)
			iz++;
	}

	// +dir_2-direction.

	for (i1 = 0; i1 < ext2; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), dir2));

		cm_eq_cm(U, U_tmp);

		if (dir2 == 0)
			it++;
		if (dir2 == 1)
			ix++;
		if (dir2 == 2)
			iy++;
		if (dir2 == 3)
			iz++;
	}

	// -dir_1-direction.

	for (i1 = 0; i1 < ext1; i1++) {
		if (dir1 == 0)
			it--;
		if (dir1 == 1)
			ix--;
		if (dir1 == 2)
			iy--;
		if (dir1 == 3)
			iz--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), dir1));

		cm_eq_cm(U, U_tmp);
	}

	// -dir_2-direction.

	for (i1 = 0; i1 < ext2; i1++) {
		if (dir2 == 0)
			it--;
		if (dir2 == 1)
			ix--;
		if (dir2 == 2)
			iy--;
		if (dir2 == 3)
			iz--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), dir2));

		cm_eq_cm(U, U_tmp);
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// gauge_field_1: spatial path at T
// gauge_field_2: spatial path at 0

complex Wilson_Loop_(double *gauge_field_1, double *gauge_field_2, int T, int L, int it, int ix, int iy, int iz, int dir2,
		int ext1, int ext2) {
	complex c1;
	int i1;

	int dir1 = 0;

	if (dir1 < 0 || dir1 > 3 || dir2 < 0 || dir2 > 3 || dir1 == dir2) {
		fprintf(stderr, "Error: complex Wilson_Loop_(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8];

	cm_eq_id(U);

	// +dir_1-direction.

	for (i1 = 0; i1 < ext1; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), dir1));

		cm_eq_cm(U, U_tmp);

		if (dir1 == 0)
			it++;
		if (dir1 == 1)
			ix++;
		if (dir1 == 2)
			iy++;
		if (dir1 == 3)
			iz++;
	}

	// +dir_2-direction.

	for (i1 = 0; i1 < ext2; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), dir2));

		cm_eq_cm(U, U_tmp);

		if (dir2 == 0)
			it++;
		if (dir2 == 1)
			ix++;
		if (dir2 == 2)
			iy++;
		if (dir2 == 3)
			iz++;
	}

	// -dir_1-direction.

	for (i1 = 0; i1 < ext1; i1++) {
		if (dir1 == 0)
			it--;
		if (dir1 == 1)
			ix--;
		if (dir1 == 2)
			iy--;
		if (dir1 == 3)
			iz--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), dir1));

		cm_eq_cm(U, U_tmp);
	}

	// -dir_2-direction.

	for (i1 = 0; i1 < ext2; i1++) {
		if (dir2 == 0)
			it--;
		if (dir2 == 1)
			ix--;
		if (dir2 == 2)
			iy--;
		if (dir2 == 3)
			iz--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field_2 + ggi(get_index(it, ix, iy, iz, T, L), dir2));

		cm_eq_cm(U, U_tmp);
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// ********************

// Computes a 2d-diagonal (off-axis) Wilson loop.

// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// (ext_T,ext_L): extension in temporal and spatial direction.

// Slow version.
complex Wilson_Loop_Diagonal_2(double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz,
		int ext_T, int ext_L) {
	complex c1;
	int i1;

	if ((dirx == 0 && ((diry == -1 || diry == +1) && (dirz == -1 || dirz == +1))
			|| diry == 0 && ((dirx == -1 || dirx == +1) && (dirz == -1 || dirz == +1))
			|| dirz == 0 && ((dirx == -1 || dirx == +1) && (diry == -1 || diry == +1))) == false) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_2(...\n");
		exit(EXIT_FAILURE);
	}

	if (ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_2(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	Diagonal_Path_2(U_tmp_, gauge_field, T, L, it, ix, iy, iz, +dirx, +diry, +dirz, ext_L);

	cm_eq_cm_ti_cm(U_tmp, U, U_tmp_);
	cm_eq_cm(U, U_tmp);

	ix += dirx * ext_L;
	iy += diry * ext_L;
	iz += dirz * ext_L;

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	Diagonal_Path_2(U_tmp_, gauge_field, T, L, it, ix, iy, iz, -dirx, -diry, -dirz, ext_L);

	cm_eq_cm_ti_cm(U_tmp, U, U_tmp_);
	cm_eq_cm(U, U_tmp);

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// Fast version.
complex Wilson_Loop_Diagonal_2_(double *gauge_field, double *diagonal_links, int T, int L, int it, int ix, int iy, int iz,
		int dirx, int diry, int dirz, int ext_T, int ext_L) {
	complex c1;
	int i1;

	if ((dirx == 0 && ((diry == -1 || diry == +1) && (dirz == -1 || dirz == +1))
			|| diry == 0 && ((dirx == -1 || dirx == +1) && (dirz == -1 || dirz == +1))
			|| dirz == 0 && ((dirx == -1 || dirx == +1) && (diry == -1 || diry == +1))) == false) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_2_(...\n");
		exit(EXIT_FAILURE);
	}

	if (ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_2_(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, diagonal_links + get_index(it, ix, iy, iz, T, L) * 8);

		cm_eq_cm(U, U_tmp);

		ix += dirx;
		iy += diry;
		iz += dirz;
	}

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		ix -= dirx;
		iy -= diry;
		iz -= dirz;

		cm_eq_cm_ti_cm_dag(U_tmp, U, diagonal_links + get_index(it, ix, iy, iz, T, L) * 8);

		cm_eq_cm(U, U_tmp);
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// ********************

// Computes a 3d-diagonal (off-axis) Wilson loop.

// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// (ext_T,ext_L): extension in temporal and spatial direction.

// Slow version.
complex Wilson_Loop_Diagonal_3(double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz,
		int ext_T, int ext_L) {
	complex c1;
	int i1;

	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != -1 && dirz != +1 || ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_3(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	Diagonal_Path_3(U_tmp_, gauge_field, T, L, it, ix, iy, iz, +dirx, +diry, +dirz, ext_L);

	cm_eq_cm_ti_cm(U_tmp, U, U_tmp_);
	cm_eq_cm(U, U_tmp);

	ix += dirx * ext_L;
	iy += diry * ext_L;
	iz += dirz * ext_L;

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	Diagonal_Path_3(U_tmp_, gauge_field, T, L, it, ix, iy, iz, -dirx, -diry, -dirz, ext_L);

	cm_eq_cm_ti_cm(U_tmp, U, U_tmp_);
	cm_eq_cm(U, U_tmp);

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// Fast version.
complex Wilson_Loop_Diagonal_3_(double *gauge_field, double *diagonal_links, int T, int L, int it, int ix, int iy, int iz,
		int dirx, int diry, int dirz, int ext_T, int ext_L) {
	complex c1;
	int i1;

	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != -1 && dirz != +1 || ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_3_(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, diagonal_links + get_index(it, ix, iy, iz, T, L) * 8);

		cm_eq_cm(U, U_tmp);

		ix += dirx;
		iy += diry;
		iz += dirz;
	}

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		ix -= dirx;
		iy -= diry;
		iz -= dirz;

		cm_eq_cm_ti_cm_dag(U_tmp, U, diagonal_links + get_index(it, ix, iy, iz, T, L) * 8);

		cm_eq_cm(U, U_tmp);
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// Replace one diagonal link by an on-axis link.
complex Wilson_Loop_Diagonal_3a_(double *gauge_field, double *diagonal_links, int T, int L, int it, int ix, int iy, int iz,
		int dirx, int diry, int dirz, int ext_T, int ext_L) {
	complex c1;
	int i1;

	// if(dirx != +1 ||
	//    diry != -1 && diry != +1 ||
	//    dirz != -1 && dirz != +1 ||
	// if(dirx != -1 && dirx != +1 ||
	//    diry != +1 ||
	//    dirz != -1 && dirz != +1 ||
	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != +1 || ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_3a_(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		if (i1 == ext_L / 2) {
			cm_eq_cm_ti_cm(U_tmp, U,
			// gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1));
			// gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2));
					gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3));

			cm_eq_cm(U, U_tmp);

			// ix++;
			// iy++;
			iz++;
		} else {
			cm_eq_cm_ti_cm(U_tmp, U, diagonal_links + get_index(it, ix, iy, iz, T, L) * 8);

			cm_eq_cm(U, U_tmp);

			ix += dirx;
			iy += diry;
			iz += dirz;
		}
	}

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		if (i1 == ext_L / 2) {
			// ix--;
			// iy--;
			iz--;

			cm_eq_cm_ti_cm_dag(U_tmp, U,
			// gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1));
			// gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2));
					gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3));

			cm_eq_cm(U, U_tmp);
		} else {
			ix -= dirx;
			iy -= diry;
			iz -= dirz;

			cm_eq_cm_ti_cm_dag(U_tmp, U, diagonal_links + get_index(it, ix, iy, iz, T, L) * 8);

			cm_eq_cm(U, U_tmp);
		}
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// Replace one diagonal link by a 2d diagonal link.
complex Wilson_Loop_Diagonal_3b_(double *gauge_field, double *diagonal_links, int T, int L, int it, int ix, int iy, int iz,
		int dirx, int diry, int dirz, int ext_T, int ext_L) {
	complex c1;
	int i1;

	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != -1 && dirz != +1 || ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_3b_(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		if (i1 == ext_L / 2) {
			double U_2d[8];
			Diagonal_Link_2(U_2d, gauge_field, T, L, it, ix, iy, iz, dirx, diry, 0);

			cm_eq_cm_ti_cm(U_tmp, U, U_2d);

			cm_eq_cm(U, U_tmp);

			ix += dirx;
			iy += diry;
		} else {
			cm_eq_cm_ti_cm(U_tmp, U, diagonal_links + get_index(it, ix, iy, iz, T, L) * 8);

			cm_eq_cm(U, U_tmp);

			ix += dirx;
			iy += diry;
			iz += dirz;
		}
	}

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		if (i1 == ext_L / 2) {
			ix -= dirx;
			iy -= diry;

			double U_2d[8];
			Diagonal_Link_2(U_2d, gauge_field, T, L, it, ix, iy, iz, dirx, diry, 0);

			cm_eq_cm_ti_cm_dag(U_tmp, U, U_2d);

			cm_eq_cm(U, U_tmp);
		} else {
			ix -= dirx;
			iy -= diry;
			iz -= dirz;

			cm_eq_cm_ti_cm_dag(U_tmp, U, diagonal_links + get_index(it, ix, iy, iz, T, L) * 8);

			cm_eq_cm(U, U_tmp);
		}
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// gauge_field_1: spatial path at T
// gauge_field_2: spatial path at 0

complex Wilson_Loop_Diagonal_3__(double *gauge_field_1, double *gauge_field_2, double *diagonal_links_1, double *diagonal_links_2,
		int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int ext_T, int ext_L) {
	complex c1;
	int i1;

	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != -1 && dirz != +1 || ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_3a__(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, diagonal_links_1 + get_index(it, ix, iy, iz, T, L) * 8);

		cm_eq_cm(U, U_tmp);

		ix += dirx;
		iy += diry;
		iz += dirz;
	}

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		ix -= dirx;
		iy -= diry;
		iz -= dirz;

		cm_eq_cm_ti_cm_dag(U_tmp, U, diagonal_links_2 + get_index(it, ix, iy, iz, T, L) * 8);

		cm_eq_cm(U, U_tmp);
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// Replace one diagonal link by an on-axis link.
complex Wilson_Loop_Diagonal_3a__(double *gauge_field_1, double *gauge_field_2, double *diagonal_links_1,
		double *diagonal_links_2, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int ext_T,
		int ext_L) {
	complex c1;
	int i1;

	// if(dirx != +1 ||
	//    diry != -1 && diry != +1 ||
	//    dirz != -1 && dirz != +1 ||
	// if(dirx != -1 && dirx != +1 ||
	//    diry != +1 ||
	//    dirz != -1 && dirz != +1 ||
	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != +1 || ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_3a__(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		if (i1 == ext_L / 2) {
			cm_eq_cm_ti_cm(U_tmp, U,
			// gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 1));
			// gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 2));
					gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 3));

			cm_eq_cm(U, U_tmp);

			// ix++;
			// iy++;
			iz++;
		} else {
			cm_eq_cm_ti_cm(U_tmp, U, diagonal_links_1 + get_index(it, ix, iy, iz, T, L) * 8);

			cm_eq_cm(U, U_tmp);

			ix += dirx;
			iy += diry;
			iz += dirz;
		}
	}

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		if (i1 == ext_L / 2) {
			// ix--;
			// iy--;
			iz--;

			cm_eq_cm_ti_cm_dag(U_tmp, U,
			// gauge_field_2 + ggi(get_index(it, ix, iy, iz, T, L), 1));
			// gauge_field_2 + ggi(get_index(it, ix, iy, iz, T, L), 2));
					gauge_field_2 + ggi(get_index(it, ix, iy, iz, T, L), 3));

			cm_eq_cm(U, U_tmp);
		} else {
			ix -= dirx;
			iy -= diry;
			iz -= dirz;

			cm_eq_cm_ti_cm_dag(U_tmp, U, diagonal_links_2 + get_index(it, ix, iy, iz, T, L) * 8);

			cm_eq_cm(U, U_tmp);
		}
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// Replace one diagonal link by a 2d diagonal link.
complex Wilson_Loop_Diagonal_3b__(double *gauge_field_1, double *gauge_field_2, double *diagonal_links_1,
		double *diagonal_links_2, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz, int ext_T,
		int ext_L) {
	complex c1;
	int i1;

	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != -1 && dirz != +1 || ext_T < 0 || ext_L < 0) {
		fprintf(stderr, "Error: complex Wilson_Loop_Diagonal_3b__(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(U);

	// positive temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		cm_eq_cm_ti_cm(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);

		it++;
	}

	// positive spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		if (i1 == ext_L / 2) {
			double U_2d[8];
			Diagonal_Link_2(U_2d, gauge_field_1, T, L, it, ix, iy, iz, dirx, diry, 0);

			cm_eq_cm_ti_cm(U_tmp, U, U_2d);

			cm_eq_cm(U, U_tmp);

			ix += dirx;
			iy += diry;
		} else {
			cm_eq_cm_ti_cm(U_tmp, U, diagonal_links_1 + get_index(it, ix, iy, iz, T, L) * 8);

			cm_eq_cm(U, U_tmp);

			ix += dirx;
			iy += diry;
			iz += dirz;
		}
	}

	// negative temporal direction

	for (i1 = 0; i1 < ext_T; i1++) {
		it--;

		cm_eq_cm_ti_cm_dag(U_tmp, U, gauge_field_1 + ggi(get_index(it, ix, iy, iz, T, L), 0));

		cm_eq_cm(U, U_tmp);
	}

	// negative spatial direction

	for (i1 = 0; i1 < ext_L; i1++) {
		if (i1 == ext_L / 2) {
			ix -= dirx;
			iy -= diry;

			double U_2d[8];
			Diagonal_Link_2(U_2d, gauge_field_2, T, L, it, ix, iy, iz, dirx, diry, 0);

			cm_eq_cm_ti_cm_dag(U_tmp, U, U_2d);

			cm_eq_cm(U, U_tmp);
		} else {
			ix -= dirx;
			iy -= diry;
			iz -= dirz;

			cm_eq_cm_ti_cm_dag(U_tmp, U, diagonal_links_2 + get_index(it, ix, iy, iz, T, L) * 8);

			cm_eq_cm(U, U_tmp);
		}
	}

	co_eq_tr_cm(&c1, U);

	c1.re /= 2.0;
	c1.im /= 2.0;

	return c1;
}

// ********************

// Computes a 1d path of spatial links.

// A: the resultung path.
// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// length: length of the path.

void Path_1(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz,
		int length) {
	if (length < 0) {
		fprintf(stderr, "Error: void Path_1(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8], U_tmp_[8];

	cm_eq_id(A);

	if (dirx == -1 || dirx == +1) {
		if (diry != 0 || dirz != 0) {
			fprintf(stderr, "Error: void Path_1(...\n");
			exit(EXIT_FAILURE);
		}

		while (1) {
			if (length == 0)
				break;

			cm_eq_cm(U_tmp, A);

			if (dirx == -1) {
				cm_eq_cm_ti_cm_dag(A, U_tmp, gauge_field + ggi(get_index(it, ix - 1, iy, iz, T, L), 1));
			} else {
				cm_eq_cm_ti_cm(A, U_tmp, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1));
			}

			ix += dirx;

			length--;
		}
	} else if (diry == -1 || diry == +1) {
		if (dirx != 0 || dirz != 0) {
			fprintf(stderr, "Error: void Path_1(...\n");
			exit(EXIT_FAILURE);
		}

		while (1) {
			if (length == 0)
				break;

			cm_eq_cm(U_tmp, A);

			if (diry == -1) {
				cm_eq_cm_ti_cm_dag(A, U_tmp, gauge_field + ggi(get_index(it, ix, iy - 1, iz, T, L), 2));
			} else {
				cm_eq_cm_ti_cm(A, U_tmp, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2));
			}

			iy += diry;

			length--;
		}
	} else if (dirz == -1 || dirz == +1) {
		if (dirx != 0 || diry != 0) {
			fprintf(stderr, "Error: void Path_1(...\n");
			exit(EXIT_FAILURE);
		}

		while (1) {
			if (length == 0)
				break;

			cm_eq_cm(U_tmp, A);

			if (dirz == -1) {
				cm_eq_cm_ti_cm_dag(A, U_tmp, gauge_field + ggi(get_index(it, ix, iy, iz - 1, T, L), 3));
			} else {
				cm_eq_cm_ti_cm(A, U_tmp, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3));
			}

			iz += dirz;

			length--;
		}
	}
}

// Computes a 1d path of temporal links.

void Path_T(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int length) {
	if (length < 0) {
		fprintf(stderr, "Error: void Path_T(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8];

	cm_eq_id(A);

	while (1) {
		if (length == 0)
			break;

		cm_eq_cm(U, A);

		cm_eq_cm_ti_cm(A, U, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0));

		it++;
		length--;
	}
}

// ********************

// Computes a 2d diagonal path of spatial links.

// A: the resultung diagonal path.
// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// length: length of the path / sqrt(2).

void Diagonal_Path_2(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz,
		int length) {
	if (length < 0) {
		fprintf(stderr, "Error: void Diagonal_Path_2(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8];

	cm_eq_id(A);

	while (1) {
		if (length == 0)
			break;

		Diagonal_Link_2(U, gauge_field, T, L, it, ix, iy, iz, dirx, diry, dirz);

		cm_eq_cm(U_tmp, A);
		cm_eq_cm_ti_cm(A, U_tmp, U);

		ix += dirx;
		iy += diry;
		iz += dirz;

		length--;
	}
}

// ********************

// Computes a 3d diagonal path of spatial links.

// A: the resultung diagonal path.
// (it,ix,iy,iz): starting point.
// (dirx,diry,dirz): spatial direction.
// length: length of the path / sqrt(3).

void Diagonal_Path_3(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz,
		int length) {
	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != -1 && dirz != +1 || length < 0) {
		fprintf(stderr, "Error: void Diagonal_Path_3(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8];

	cm_eq_id(A);

	while (1) {
		if (length == 0)
			break;

		Diagonal_Link_3(U, gauge_field, T, L, it, ix, iy, iz, dirx, diry, dirz);

		cm_eq_cm(U_tmp, A);
		cm_eq_cm_ti_cm(A, U_tmp, U);

		ix += dirx;
		iy += diry;
		iz += dirz;

		length--;
	}
}

// Replace one diagonal link by an on-axis link.
void Diagonal_Path_3a(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz,
		int length) {
	int i1;

	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != +1 || length < 0) {
		fprintf(stderr, "Error: void Diagonal_Path_3a(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8];

	cm_eq_id(A);

	for (i1 = 0; i1 < length; i1++) {
		if (i1 == length / 2) {
			cm_eq_cm(U_tmp, A);
			cm_eq_cm_ti_cm(A, U_tmp, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3));

			cm_eq_cm(U, U_tmp);

			iz++;
		} else {
			Diagonal_Link_3(U, gauge_field, T, L, it, ix, iy, iz, dirx, diry, dirz);

			cm_eq_cm(U_tmp, A);
			cm_eq_cm_ti_cm(A, U_tmp, U);

			ix += dirx;
			iy += diry;
			iz += dirz;
		}
	}
}

// Replace one diagonal link by a 2d diagonal link.
void Diagonal_Path_3b(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz,
		int length) {
	int i1;

	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != -1 && dirz != +1 || length < 0) {
		fprintf(stderr, "Error: void Diagonal_Path_3b(...\n");
		exit(EXIT_FAILURE);
	}

	double U[8], U_tmp[8];

	cm_eq_id(A);

	for (i1 = 0; i1 < length; i1++) {
		if (i1 == length / 2) {
			Diagonal_Link_2(U, gauge_field, T, L, it, ix, iy, iz, dirx, diry, 0);

			cm_eq_cm(U_tmp, A);
			cm_eq_cm_ti_cm(A, U_tmp, U);

			cm_eq_cm(U, U_tmp);

			ix += dirx;
			iy += diry;
		} else {
			Diagonal_Link_3(U, gauge_field, T, L, it, ix, iy, iz, dirx, diry, dirz);

			cm_eq_cm(U_tmp, A);
			cm_eq_cm_ti_cm(A, U_tmp, U);

			ix += dirx;
			iy += diry;
			iz += dirz;
		}
	}
}

// ********************

// Computes a 2d diagonal spatial link.

void Diagonal_Link_2(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz) {
	double U_tmp[8], U_tmp_[8];

	if (dirx == 0) {
		if (diry != -1 && diry != +1 || dirz != -1 && dirz != +1) {
			fprintf(stderr, "Error: void Diagonal_Link_2(...\n");
			exit(EXIT_FAILURE);
		}

		cm_eq_zero(A);

		// y --> z

		cm_eq_id(U_tmp_);

		if (diry == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy - 1, iz, T, L), 2));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2));
		}

		cm_eq_cm(U_tmp_, U_tmp);

		if (dirz == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy + diry, iz - 1, T, L), 3));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy + diry, iz, T, L), 3));
		}

		cm_pl_eq_cm(A, U_tmp);

		// z --> y

		cm_eq_id(U_tmp_);

		if (dirz == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz - 1, T, L), 3));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3));
		}

		cm_eq_cm(U_tmp_, U_tmp);

		if (diry == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy - 1, iz + dirz, T, L), 2));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz + dirz, T, L), 2));
		}

		cm_pl_eq_cm(A, U_tmp);

		cm_proj(A);
	} else if (diry == 0) {
		if (dirx != -1 && dirx != +1 || dirz != -1 && dirz != +1) {
			fprintf(stderr, "Error: void Diagonal_Link_2(...\n");
			exit(EXIT_FAILURE);
		}

		cm_eq_zero(A);

		// x --> z

		cm_eq_id(U_tmp_);

		if (dirx == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy, iz, T, L), 1));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1));
		}

		cm_eq_cm(U_tmp_, U_tmp);

		if (dirz == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy, iz - 1, T, L), 3));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy, iz, T, L), 3));
		}

		cm_pl_eq_cm(A, U_tmp);

		// z --> x

		cm_eq_id(U_tmp_);

		if (dirz == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz - 1, T, L), 3));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3));
		}

		cm_eq_cm(U_tmp_, U_tmp);

		if (dirx == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy, iz + dirz, T, L), 1));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz + dirz, T, L), 1));
		}

		cm_pl_eq_cm(A, U_tmp);

		cm_proj(A);
	} else if (dirz == 0) {
		if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1) {
			fprintf(stderr, "Error: void Diagonal_Link_2(...\n");
			exit(EXIT_FAILURE);
		}

		cm_eq_zero(A);

		// x --> y

		cm_eq_id(U_tmp_);

		if (dirx == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy, iz, T, L), 1));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1));
		}

		cm_eq_cm(U_tmp_, U_tmp);

		if (diry == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy - 1, iz, T, L), 2));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy, iz, T, L), 2));
		}

		cm_pl_eq_cm(A, U_tmp);

		// y --> x

		cm_eq_id(U_tmp_);

		if (diry == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy - 1, iz, T, L), 2));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2));
		}

		cm_eq_cm(U_tmp_, U_tmp);

		if (dirx == -1) {
			cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy + diry, iz, T, L), 1));
		} else {
			cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy + diry, iz, T, L), 1));
		}

		cm_pl_eq_cm(A, U_tmp);

		cm_proj(A);
	} else {
		fprintf(stderr, "Error: void Diagonal_Link_2(...\n");
		exit(EXIT_FAILURE);
	}
}

// Computes 2d diagonal spatial links.

void Diagonal_Links_2(double *diagonal_links, double *gauge_field, int T, int L, int dirx, int diry, int dirz) {
	int it, ix, iy, iz;
	double U[8];

	for (it = 0; it < T; it++) {
		for (ix = 0; ix < L; ix++) {
			for (iy = 0; iy < L; iy++) {
				for (iz = 0; iz < L; iz++) {
					Diagonal_Link_2(U, gauge_field, T, L, it, ix, iy, iz, dirx, diry, dirz);

					cm_eq_cm(diagonal_links + get_index(it, ix, iy, iz, T, L) * 8, U);
				}
			}
		}
	}
}

// ********************

// Computes a 3d diagonal spatial link.

void Diagonal_Link_3(double *A, double *gauge_field, int T, int L, int it, int ix, int iy, int iz, int dirx, int diry, int dirz) {
	double U_tmp[8], U_tmp_[8];

	if (dirx != -1 && dirx != +1 || diry != -1 && diry != +1 || dirz != -1 && dirz != +1) {
		fprintf(stderr, "Error: void Diagonal_Link_3(...\n");
		exit(EXIT_FAILURE);
	}

	cm_eq_zero(A);

	// x --> y --> z

	cm_eq_id(U_tmp_);

	if (dirx == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy, iz, T, L), 1));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (diry == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy - 1, iz, T, L), 2));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy, iz, T, L), 2));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (dirz == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy + diry, iz - 1, T, L), 3));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy + diry, iz, T, L), 3));
	}

	cm_pl_eq_cm(A, U_tmp);

	// x --> z --> y

	cm_eq_id(U_tmp_);

	if (dirx == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy, iz, T, L), 1));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 1));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (dirz == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy, iz - 1, T, L), 3));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy, iz, T, L), 3));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (diry == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy - 1, iz + dirz, T, L), 2));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy, iz + dirz, T, L), 2));
	}

	cm_pl_eq_cm(A, U_tmp);

	// y --> x --> z

	cm_eq_id(U_tmp_);

	if (diry == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy - 1, iz, T, L), 2));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (dirx == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy + diry, iz, T, L), 1));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy + diry, iz, T, L), 1));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (dirz == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy + diry, iz - 1, T, L), 3));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy + diry, iz, T, L), 3));
	}

	cm_pl_eq_cm(A, U_tmp);

	// y --> z --> x

	cm_eq_id(U_tmp_);

	if (diry == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy - 1, iz, T, L), 2));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 2));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (dirz == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy + diry, iz - 1, T, L), 3));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy + diry, iz, T, L), 3));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (dirx == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy + diry, iz + dirz, T, L), 1));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy + diry, iz + dirz, T, L), 1));
	}

	cm_pl_eq_cm(A, U_tmp);

	// z --> x --> y

	cm_eq_id(U_tmp_);

	if (dirz == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz - 1, T, L), 3));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (dirx == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy, iz + dirz, T, L), 1));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz + dirz, T, L), 1));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (diry == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy - 1, iz + dirz, T, L), 2));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix + dirx, iy, iz + dirz, T, L), 2));
	}

	cm_pl_eq_cm(A, U_tmp);

	// z --> y --> x

	cm_eq_id(U_tmp_);

	if (dirz == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz - 1, T, L), 3));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 3));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (diry == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy - 1, iz + dirz, T, L), 2));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy, iz + dirz, T, L), 2));
	}

	cm_eq_cm(U_tmp_, U_tmp);

	if (dirx == -1) {
		cm_eq_cm_ti_cm_dag(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix - 1, iy + diry, iz + dirz, T, L), 1));
	} else {
		cm_eq_cm_ti_cm(U_tmp, U_tmp_, gauge_field + ggi(get_index(it, ix, iy + diry, iz + dirz, T, L), 1));
	}

	cm_pl_eq_cm(A, U_tmp);

	cm_proj(A);
}

// Computes 3d diagonal spatial links.

void Diagonal_Links_3(double *diagonal_links, double *gauge_field, int T, int L, int dirx, int diry, int dirz) {
	int it, ix, iy, iz;
	double U[8];

	for (it = 0; it < T; it++) {
		for (ix = 0; ix < L; ix++) {
			for (iy = 0; iy < L; iy++) {
				for (iz = 0; iz < L; iz++) {
					Diagonal_Link_3(U, gauge_field, T, L, it, ix, iy, iz, dirx, diry, dirz);

					cm_eq_cm(diagonal_links + get_index(it, ix, iy, iz, T, L) * 8, U);
				}
			}
		}
	}
}

// ********************

// Computes the average plaquette.

double Average_Plaquette(double *gauge_field, int T, int L) {
	int i1, i2;

	int it, ix, iy, iz;

	double P = 0.0;

	for (it = 0; it < T; it++) {
		for (ix = 0; ix < L; ix++) {
			for (iy = 0; iy < L; iy++) {
				for (iz = 0; iz < L; iz++) {
					for (i1 = 0; i1 < 3; i1++) {
						for (i2 = i1 + 1; i2 < 4; i2++) {
							complex c1 = Wilson_Loop(gauge_field, T, L, it, ix, iy, iz, i1, i2, 1, 1);

							P += c1.re;
						}
					}
				}
			}
		}
	}

	return P / (double) (T * L * L * L * 6);
}

// ********************
