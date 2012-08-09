#include <gsl/gsl_vector.h>
#include "basis.h"
#include "common.h"

#ifndef __INTEGRAL__INT1E__
#define __INTEGRAL__INT1E__
double RecCoeff(int i, int j, int t, double *alphpP, double *PA, double *PB);
double R(int n, int t, int u, int v, const gsl_vector *PX, double gamma, FMW *, int );
double overlap_gto(const GTO* g1, const gsl_vector* A, const GTO* g2, const gsl_vector* B, int debug);
double nuclear_elect_attraction_gto(const GTO* g1, const gsl_vector* A, \
        const GTO* g2, const gsl_vector* B, const gsl_vector* C, FMW *fmw, int debug);
#endif
