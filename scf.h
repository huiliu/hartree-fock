#include <gsl/gsl_matrix.h>

#ifndef __INTEGRAL__SCF__
#define __INTEGRAL__SCF__

void HartreeFock();
gsl_matrix* S_i_root(gsl_matrix *S, int n);
gsl_matrix *Fock(gsl_matrix *h, double ****e2_int, double** density, int n);
void density(double **, gsl_matrix* coef, int n);
gsl_matrix* scf(gsl_matrix *f, const gsl_matrix *s_root, int n, double* energy);

#endif
