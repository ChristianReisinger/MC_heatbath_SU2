// geometry.hh

#ifndef __GEOMETRY_HH__

#define __GEOMETRY_HH__

// ********************

inline int get_index(const int t, const int x, const int y, const int z, const int T, const int L) {
	int tt = (t + T) % T;
	int xx = (x + L) % L;
	int yy = (y + L) % L;
	int zz = (z + L) % L;

	return ((tt * L + xx) * L + yy) * L + zz;
}

inline int get_index_timeslice(const int x, const int y, const int z, const int T, const int L) {
	int xx = (x + L) % L;
	int yy = (y + L) % L;
	int zz = (z + L) % L;

	return (xx * L + yy) * L + zz;
}

inline int get_index_timeslice_t(const int t, const int x, const int y, const int z, const int T, const int L) {
	int xx = (x + L) % L;
	int yy = (y + L) % L;
	int zz = (z + L) % L;

	return (xx * L + yy) * L + zz;
}

inline int ggi(const int ix, const int mu) {
	return (4 * ix + mu) * 8;
}

inline int ggi2(const int ix, const int mu) {
	return (3 * ix + mu);
}

inline int get_index_WL(const int pla, const int con, const int WL_R) {
	return (WL_R * con + pla);
}

inline int get_index_plaquette2D(const int con, const int pla_parallel, const int pla_vertical1, const int WL_R,
		const int R_parallel, const int R_vertical) {
	return (con * (WL_R + 2 * R_parallel) + (pla_parallel + R_parallel)) * (R_vertical * 2 + 1) + (pla_vertical1 + R_vertical);
}

inline int get_index_res(const int pla_parallel, const int pla_vertical1, const int WL_R, const int R_parallel,
		const int R_vertical) {
	return ((pla_vertical1) * ((WL_R + 1) / 2 + R_parallel) + (pla_parallel + R_parallel));
}

inline int get_index_side(const int con, const int pla_parallel, const int WL_R, const int R_parallel) {
	return (con * (R_parallel + (WL_R - 1) / 2) + (pla_parallel + R_parallel));
}

inline int get_index_updown(const int con, const int pla_vertical1, const int R_vertical) {
	return (con * R_vertical + (pla_vertical1 - 1));
}

inline int get_index_sideupdown(const int con, const int pla_parallel, const int pla_vertical1, const int WL_R,
		const int R_parallel, const int R_vertical) {
	return ((con * (R_parallel + (WL_R - 1) / 2) + (pla_parallel + R_parallel)) * R_vertical + (pla_vertical1 - 1));
}

#endif
