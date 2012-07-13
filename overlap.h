#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "basis.h"

#ifndef __INTEGRAL__OVERLAP__
#define __INTEGRAL__OVERLAP__

int factorial(int);
int factorial_2(int);
double fact_l_lambda(int , int );
double fi_l_ll_pax_pbx(int ii, int l, int ll, double pax, double pbx, int);
double I_xyz(int l, double pax, int ll, double pbx, double gamma, int flags);
double gauss_K(double a, const gsl_vector *A, double b, const gsl_vector *B);
gsl_vector* gaussian_product_center(const double a, const gsl_vector *A, const double b, const gsl_vector *B, int);
double overlap_basis(const BASIS *, const gsl_vector *, const BASIS *, const gsl_vector *, int);
double overlap_gauss(const GTO g1, const gsl_vector* A, const GTO g2, const gsl_vector* B, int debug);
#endif
