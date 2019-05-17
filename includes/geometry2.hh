#include <vector>
#include "geometry.hh"
#include "global_defs.hh"

#ifndef _GEOMETRY2_HH_
#define _GEOMETRY2_HH_

/**
 * @param n lattice site {t,x,y,z}
 * @param T
 * @param L
 */
inline unsigned long int get_index(const std::vector<int>& n, const int T, const int L) {
	return get_index(n[0], n[1], n[2], n[3], T, L);
}

/**
 * @param p lattice position, given by the array {t,x,y,z}.
 * @param n number of SU(N) matrices on each lattice site.
 * @param i index to obtain the ith matrix on the lattice site p.
 * @return superindex of SU(N) matrices on the lattice with n matrices on each lattice site.
 */
inline unsigned long int ggi_n(const int (&p)[4], const int n, const int i, const int T, const int L) {
	return (n * get_index(p[0], p[1], p[2], p[3], T, L) + i) * (unsigned long int) SUN_elems;
}

/**
 * @param t
 * @param x
 * @param y
 * @param z
 * @param n number of SU(N) matrices on each lattice site.
 * @param i index to obtain the ith matrix on the lattice site p.
 * @return superindex of SU(N) matrices on the lattice with n matrices on each lattice site.
 */
inline unsigned long int ggi_n(const int t, const int x, const int y, const int z, const int n, const int i, const int T,
		const int L) {
	return (n * get_index(t, x, y, z, T, L) + i) * (unsigned long int) SUN_elems;
}

#endif
