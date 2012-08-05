#include "common.h"
#include "basis.h"

#ifndef __INTEGRAL__HAMILTONIAN__
#define __INTEGRAL__HAMILTONIAN__
double kinetic_I_xyz(const GTO* g1, const gsl_vector* A, 
                     const GTO* g2, const gsl_vector* B, int flags, int debug);
double kinetic_I_xyz_c(const GTO* g1, const gsl_vector* A, 
                     const GTO* g2, const gsl_vector* B, int flags, int debug);
double kinetic_gto(const GTO *g1, const gsl_vector* A, 
                   const GTO* g2, const gsl_vector* B, int debug);
double check_kinetic(const BASIS* b1, const BASIS* b2, int debug);
double kinetic_basis(const BASIS* b1, const BASIS* b2, int debug);
gsl_matrix* kinetic_matrix(INPUT_INFO* b);


double check_nuclear(const BASIS* b1, const BASIS* b2, ATOM_INFO **atomList, int atomCount);
double* A_iru(int l1, int l2, double Ax, double Bx, double Cx, double gamma, int debug);
//double nuclear_elect_attraction_gto(const GTO* g1, const gsl_vector* A, const GTO* g2, const gsl_vector* B, const gsl_vector* C, int debug);


double nuclear_elect_attraction_basis(const BASIS* b1, const BASIS* b2, 
                                ATOM_INFO **atomList, int atomCount,F_INC_GAMMA *, int debug);
gsl_matrix* nuclear_attraction_matrix(INPUT_INFO* b);
gsl_matrix* hamiltonian(INPUT_INFO* b);
#endif
