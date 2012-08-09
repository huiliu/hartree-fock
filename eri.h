#include "common.h"
#include "basis.h"

#ifndef __INTEGRAL__ERI_INCLUDE__
#define __INTEGRAL__ERI_INCLUDE__
double ERI_gto(const GTO* g1, const gsl_vector* A, 
                              const GTO* g2, const gsl_vector* B,
                              const GTO* g3, const gsl_vector* C,
                              const GTO* g4, const gsl_vector* D, 
                              FMW *, int debug);
#endif
