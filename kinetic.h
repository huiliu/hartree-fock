
#include "common.h"
#include "basis.h"

#ifndef __INTEGRAL__KINETIC__
#define __INTEGRAL__KINETIC__
int gtoIsNeg(const GTO* g);
double kinetic_I_xyz(const GTO* g1, const gsl_vector* A, 
                     const GTO* g2, gsl_vector* B, int flags, int debug);
double kinetic_I_xyz_c(const GTO* g1, const gsl_vector* A, 
                     const GTO* g2, gsl_vector* B, int flags, int debug);
double kinetic_gto(const GTO *g1, const gsl_vector* A, 
                   const GTO* g2, gsl_vector* B, int debug);
double check_kinetic(const BASIS* b1, const BASIS* b2, int debug);
double kinetic_basis(const BASIS* b1, const BASIS* b2, int debug);
gsl_matrix* kinetic_matrix(INPUT_INFO* b);
#endif
