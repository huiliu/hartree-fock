#include <gsl/gsl_vector.h>
#include "basis.h"
#include "common.h"

//double RecCoeff(int i, int j, int t, const double *alphpP, const double *PA, const double *PB);
double RecCoeff(int i, int j, int t, double *alphpP, double *PA, double *PB);
double overlap_gto(const GTO* g1, const gsl_vector* A, const GTO* g2, const gsl_vector* B, int debug);
