// linear_algebra.hh

#ifndef __LINEAR_ALGEBRA_HH__

#define __LINEAR_ALGEBRA_HH__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ********************

// A complex number.

class complex {
public:

	double re;
	double im;
};

inline void co_ti_eq_co(complex *c1, const complex *c2) {
	double c1_re = c1->re * c2->re - c1->im * c2->im;
	c1->im = c1->re * c2->im + c1->im * c2->re;
	c1->re = c1_re;
}

inline void co_eq_zero(complex *c) {
	c->re = 0.0;
	c->im = 0.0;
}

inline void co_pl_eq_co(complex *c1, const complex *c2) {
	c1->re += c2->re;
	c1->im += c2->im;
}

inline void co_di_eq_re(complex *c, double d) {
	c->re /= d;
	c->im /= d;
}

inline void co_eq_co_ti_co(complex *c1, const complex *c2, const complex *c3) {
	c1->re = c2->re * c3->re - c2->im * c3->im;
	c1->im = c2->re * c3->im + c2->im * c3->re;
}

// ********************

// A = 0.
inline void cm_eq_zero(double *A);

// A = 1.
inline void cm_eq_id(double *A);

// A = B.
inline void cm_eq_cm(double *A, const double *B);

// A = A^\dagger.
inline void cm_dag_eq_cm(double *A);

// A = B^\dagger.
inline void cm_eq_cm_dag(double *A, const double *B);

// c = det(A).
inline void co_eq_det_cm(complex *c, const double *A);

// c = tr(A).
inline void co_eq_tr_cm(complex *c, const double *A);

// A = Proj_SU2(A).
inline void cm_proj(double *A);

// Prints A.
inline void cm_fprintf(const double *A, FILE *file);

// A = A + B.
inline void cm_pl_eq_cm(double *A, const double *B);

// A = A + B^\dagger.
inline void cm_pl_eq_cm_dag(double *A, const double *B);

// A = B * d.
inline void cm_eq_cm_ti_re(double *A, const double *B, double d);

// A = A * d.
inline void cm_ti_eq_re(double *A, double d);

// A = B * c.
inline void cm_eq_cm_ti_co(double *A, const double *B, const complex *c);

// A = B * C.
inline void cm_eq_cm_ti_cm(double *A, const double *B, const double *C);

// A = B^\dagger * C.
inline void cm_eq_cm_dag_ti_cm(double *A, const double *B, const double *C);

// A = B * C^\dagger.
inline void cm_eq_cm_ti_cm_dag(double *A, const double *B, const double *C);

// A = B^\dagger * C^\dagger.
inline void cm_eq_cm_dag_ti_cm_dag(double *A, const double *B, const double *C);

// A = h[0] + i h[1] sigma[1] + i h[2] sigma[2] + i h[3] sigma[3]. 
inline void cm_from_h(double *A, const double *h);
inline void h_from_cm(double *h, const double *A);

// ********************

// A = 0.

inline void cm_eq_zero(double *A) {
	A[0] = 0.0;
	A[1] = 0.0;
	A[2] = 0.0;
	A[3] = 0.0;

	A[4] = 0.0;
	A[5] = 0.0;
	A[6] = 0.0;
	A[7] = 0.0;
}

// ********************

// A = 1.

inline void cm_eq_id(double *A) {
	A[0] = 1.0;
	A[1] = 0.0;
	A[2] = 0.0;
	A[3] = 0.0;

	A[4] = 0.0;
	A[5] = 0.0;
	A[6] = 1.0;
	A[7] = 0.0;
}

// ********************

// A = B.

inline void cm_eq_cm(double *A, const double *B) {
	A[0] = B[0];
	A[1] = B[1];
	A[2] = B[2];
	A[3] = B[3];

	A[4] = B[4];
	A[5] = B[5];
	A[6] = B[6];
	A[7] = B[7];
}

// ********************

// A = A^\dagger.

inline void cm_dag_eq_cm(double *A) {
	double d1;

	A[1] = -A[1];

	d1 = A[4];
	A[4] = A[2];
	A[2] = d1;

	d1 = -A[5];
	A[5] = -A[3];
	A[3] = d1;

	A[7] = -A[7];
}

// ********************

// A = B^\dagger.

inline void cm_eq_cm_dag(double *A, const double *B) {
	A[0] = B[0];
	A[1] = -B[1];
	A[2] = B[4];
	A[3] = -B[5];

	A[4] = B[2];
	A[5] = -B[3];
	A[6] = B[6];
	A[7] = -B[7];
}

// ********************

// c = det(A).

inline void co_eq_det_cm(complex *c, const double *A) {
	c->re = A[0] * A[6] - A[1] * A[7] - A[2] * A[4] + A[3] * A[5];
	c->im = A[0] * A[7] + A[1] * A[6] - A[2] * A[5] - A[3] * A[4];
}

// ********************

// c = tr(A).

inline void co_eq_tr_cm(complex *c, const double *A) {
	c->re = A[0] + A[6];
	c->im = A[1] + A[7];
}

// ********************

// A = Proj_SU2(A).
// Projects a matrix of the form A = (a, b; -b*, a*) on SU(2).
// P_{SU(2)}(A) = A / sqrt(det(A)).
inline void cm_proj(double *A) {
	double _1_sqrt_det =
			1.0 / sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2] + A[3] * A[3]);

	cm_ti_eq_re(A, _1_sqrt_det);

	/*
	 // !!!!!!!!!!
	 // !!!!!!!!!!
	 // !!!!!!!!!!
	 // This is too complicated!
	 // !!!!!!!!!!
	 // !!!!!!!!!!
	 // !!!!!!!!!!


	 // Compute the determinant.

	 complex det;
	 co_eq_det_cm(&det, A);


	 // Compute 1 / sqrt(det) (take the root closest to 1.0).

	 double _1_sqrt_r = 1.0 / pow(det.re*det.re + det.im*det.im, 0.25);
	 double phi_2 = 0.5 * atan2(det.im, det.re);

	 complex _1_sqrt_det;
	 _1_sqrt_det.re = _1_sqrt_r * cos(+phi_2);
	 _1_sqrt_det.im = _1_sqrt_r * sin(-phi_2);


	 // Project the matrix.

	 double A_proj[8];
	 cm_eq_cm_ti_co(A_proj, A, &_1_sqrt_det);


	 // Copy the result.

	 cm_eq_cm(A, A_proj);
	 */
}

// ********************

// Prints A.

inline void cm_fprintf(const double *A, FILE *file) {
	fprintf(file, "%+9.6lf %+9.6lf I   %+9.6lf %+9.6lf I\n",
			A[0], A[1], A[2], A[3]);

	fprintf(file, "%+9.6lf %+9.6lf I   %+9.6lf %+9.6lf I\n",
			A[4], A[5], A[6], A[7]);
}

// ********************

// A = A + B.

inline void cm_pl_eq_cm(double *A, const double *B) {
	A[0] += B[0];
	A[1] += B[1];
	A[2] += B[2];
	A[3] += B[3];

	A[4] += B[4];
	A[5] += B[5];
	A[6] += B[6];
	A[7] += B[7];
}

// ********************

// A = A + B^\dagger.

inline void cm_pl_eq_cm_dag(double *A, const double *B) {
	A[0] += B[0];
	A[1] -= B[1];
	A[2] += B[4];
	A[3] -= B[5];

	A[4] += B[2];
	A[5] -= B[3];
	A[6] += B[6];
	A[7] -= B[7];
}

// ********************

// A = B * d.

inline void cm_eq_cm_ti_re(double *A, const double *B, double d) {
	A[0] = B[0] * d;
	A[1] = B[1] * d;
	A[2] = B[2] * d;
	A[3] = B[3] * d;

	A[4] = B[4] * d;
	A[5] = B[5] * d;
	A[6] = B[6] * d;
	A[7] = B[7] * d;
}

// ********************

// A = A * d.

inline void cm_ti_eq_re(double *A, double d) {
	A[0] *= d;
	A[1] *= d;
	A[2] *= d;
	A[3] *= d;

	A[4] *= d;
	A[5] *= d;
	A[6] *= d;
	A[7] *= d;
}

// ********************

// A = B * c.

inline void cm_eq_cm_ti_co(double *A, const double *B, const complex *c) {
	double re = c->re;
	double im = c->im;

	A[0] = B[0] * re - B[1] * im;
	A[1] = B[1] * re + B[0] * im;
	A[2] = B[2] * re - B[3] * im;
	A[3] = B[3] * re + B[2] * im;

	A[4] = B[4] * re - B[5] * im;
	A[5] = B[5] * re + B[4] * im;
	A[6] = B[6] * re - B[7] * im;
	A[7] = B[7] * re + B[6] * im;
}

// ********************

// A = B * C.

inline void cm_eq_cm_ti_cm(double *A, const double *B, const double *C) {
	// A00 = B00*C00 + B01*C10
	A[0] = B[0] * C[0] - B[1] * C[1] + B[2] * C[4] - B[3] * C[5];
	A[1] = B[0] * C[1] + B[1] * C[0] + B[2] * C[5] + B[3] * C[4];

	// A01 = B00*C01 + B01*C11
	A[2] = B[0] * C[2] - B[1] * C[3] + B[2] * C[6] - B[3] * C[7];
	A[3] = B[0] * C[3] + B[1] * C[2] + B[2] * C[7] + B[3] * C[6];

	// A10 = B10*C00 + B11*C10
	A[4] = B[4] * C[0] - B[5] * C[1] + B[6] * C[4] - B[7] * C[5];
	A[5] = B[4] * C[1] + B[5] * C[0] + B[6] * C[5] + B[7] * C[4];

	// A11 = B10*C01 + B11*C11
	A[6] = B[4] * C[2] - B[5] * C[3] + B[6] * C[6] - B[7] * C[7];
	A[7] = B[4] * C[3] + B[5] * C[2] + B[6] * C[7] + B[7] * C[6];
}

// ********************

// A = B^\dagger * C.

inline void cm_eq_cm_dag_ti_cm(double *A, const double *B, const double *C) {
	// A00 = B00*C00 + B01*C10
	A[0] =
			B[0] * C[0] + B[1] * C[1] +
					B[4] * C[4] + B[5] * C[5];
	A[1] =
			B[0] * C[1] - B[1] * C[0] +
					B[4] * C[5] - B[5] * C[4];

	// A01 = B00*C01 + B01*C11
	A[2] =
			B[0] * C[2] + B[1] * C[3] +
					B[4] * C[6] + B[5] * C[7];
	A[3] =
			B[0] * C[3] - B[1] * C[2] +
					B[4] * C[7] - B[5] * C[6];

	// A10 = B10*C00 + B11*C10
	A[4] =
			B[2] * C[0] + B[3] * C[1] +
					B[6] * C[4] + B[7] * C[5];
	A[5] =
			B[2] * C[1] - B[3] * C[0] +
					B[6] * C[5] - B[7] * C[4];

	// A11 = B10*C01 + B11*C11
	A[6] =
			B[2] * C[2] + B[3] * C[3] +
					B[6] * C[6] + B[7] * C[7];
	A[7] =
			B[2] * C[3] - B[3] * C[2] +
					B[6] * C[7] - B[7] * C[6];
}

// ********************

// A = B * C^\dagger.

inline void cm_eq_cm_ti_cm_dag(double *A, const double *B, const double *C) {
	// A00 = B00*C00 + B01*C10
	A[0] =
			B[0] * C[0] + B[1] * C[1] +
					B[2] * C[2] + B[3] * C[3];
	A[1] =
			-B[0] * C[1] + B[1] * C[0] +
					-B[2] * C[3] + B[3] * C[2];

	// A01 = B00*C01 + B01*C11
	A[2] =
			B[0] * C[4] + B[1] * C[5] +
					B[2] * C[6] + B[3] * C[7];
	A[3] =
			-B[0] * C[5] + B[1] * C[4] +
					-B[2] * C[7] + B[3] * C[6];

	// A10 = B10*C00 + B11*C10
	A[4] =
			B[4] * C[0] + B[5] * C[1] +
					B[6] * C[2] + B[7] * C[3];
	A[5] =
			-B[4] * C[1] + B[5] * C[0] +
					-B[6] * C[3] + B[7] * C[2];

	// A11 = B10*C01 + B11*C11
	A[6] =
			B[4] * C[4] + B[5] * C[5] +
					B[6] * C[6] + B[7] * C[7];
	A[7] =
			-B[4] * C[5] + B[5] * C[4] +
					-B[6] * C[7] + B[7] * C[6];
}

// ********************

// A = B^\dagger * C^\dagger.

inline void cm_eq_cm_dag_ti_cm_dag(double *A, const double *B, const double *C) {
	// A00 = B00*C00 + B01*C10
	A[0] =
			B[0] * C[0] - B[1] * C[1] +
					B[4] * C[2] - B[5] * C[3];
	A[1] =
			-B[0] * C[1] - B[1] * C[0] +
					-B[4] * C[3] - B[5] * C[2];

	// A01 = B00*C01 + B01*C11
	A[2] =
			B[0] * C[4] - B[1] * C[5] +
					B[4] * C[6] - B[5] * C[7];
	A[3] =
			-B[0] * C[5] - B[1] * C[4] +
					-B[4] * C[7] - B[5] * C[6];

	// A10 = B10*C00 + B11*C10
	A[4] =
			B[2] * C[0] - B[3] * C[1] +
					B[6] * C[2] - B[7] * C[3];
	A[5] =
			-B[2] * C[1] - B[3] * C[0] +
					-B[6] * C[3] - B[7] * C[2];

	// A11 = B10*C01 + B11*C11
	A[6] =
			B[2] * C[4] - B[3] * C[5] +
					B[6] * C[6] - B[7] * C[7];
	A[7] =
			-B[2] * C[5] - B[3] * C[4] +
					-B[6] * C[7] - B[7] * C[6];
}

// ********************

// A = h[0] + i h[1] sigma[1] + i h[2] sigma[2] + i h[3] sigma[3]. 

inline void cm_from_h(double *A, const double *h) {
	A[0] = +h[0];
	A[1] = +h[3];
	A[2] = +h[2];
	A[3] = +h[1];

	A[4] = -h[2];
	A[5] = +h[1];
	A[6] = +h[0];
	A[7] = -h[3];
}

inline void h_from_cm(double *h, const double *A) {
	h[0] = A[0];
	h[3] = A[1];
	h[2] = A[2];
	h[1] = A[3];
}

#endif
