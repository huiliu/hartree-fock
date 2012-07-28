#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include "basis.h"

#ifndef __INTEGRAL__OVERLAP__
#define __INTEGRAL__OVERLAP__

double fact_l_lambda(int , int );
double fi_l_ll_pax_pbx(int ii, int l, int ll, double pax, double pbx, int);
double I_xyz(int l, double pax, int ll, double pbx, double gamma, int flags);
double gauss_K(double a, const gsl_vector *A, double b, const gsl_vector *B);
gsl_vector* gaussian_product_center(const double a, const gsl_vector *A, const double b, const gsl_vector *B, int);
double overlap_basis(const BASIS *, const gsl_vector *, const BASIS *, const gsl_vector *, int);
double overlap_gto(const GTO*, const gsl_vector* A, const GTO*, const gsl_vector* B, int debug);
gsl_matrix* overlap_matrix(INPUT_INFO* );
#endif
