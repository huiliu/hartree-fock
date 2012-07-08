#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "basis.h"

int factorial(int);
int factorial_2(int);
double fact_l_lambda(int , int );
double fi_l_ll_pax_pbx(int ii, int l, int ll, double pax, double pbx);
double I_xyz(double a, int l, double pax, double b, int ll, double pbx);
double gauss_K(double a, const gsl_vector *A, double b, const gsl_vector *B);
gsl_vector* gaussian_product_center(const double a, const gsl_vector *A, const double b, const gsl_vector *B);
double normalize_coeff(const GTO *g);
double overlap_single(const GTO *g1, const gsl_vector *A, const GTO *g2, const gsl_vector *B);
