#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basis.h"
#include "overlap.h"
#include "common.h"
#include "ints.h"

#ifndef __INTEGRAL__NUCLEAR__
#define __INTEGRAL__NUCLEAR__
double check_nuclear(const BASIS* b1, const BASIS* b2, ATOM_INFO **atomList, int atomCount);
double* A_iru(int l1, int l2, double Ax, double Bx, double Cx, double gamma, int debug);
double nuclear_elect_attraction_gto(const GTO* g1, const gsl_vector* A, \
        const GTO* g2, const gsl_vector* B, const gsl_vector* C, int debug);


double nuclear_elect_attraction_basis(const BASIS* b1, const BASIS* b2, 
                                ATOM_INFO **atomList, int atomCount, int debug);
gsl_matrix* nuclear_attraction_matrix(INPUT_INFO* b);
#endif
