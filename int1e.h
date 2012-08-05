#include <gsl/gsl_vector.h>
#include "basis.h"
#include "common.h"

//double RecCoeff(int i, int j, int t, const double *alphpP, const double *PA, const double *PB);
double RecCoeff(int i, int j, int t, double *alphpP, double *PA, double *PB);
//double R(int n, int t, int u, int v, const gsl_vector *PX, double gamma, int debug);
double R(int n, int t, int u, int v, const gsl_vector *PX, double gamma, F_INC_GAMMA *f, int debug);
double overlap_gto(const GTO* g1, const gsl_vector* A, const GTO* g2, const gsl_vector* B, int debug);
double nuclear_elect_attraction_gto(const GTO* g1, const gsl_vector* A, \
        const GTO* g2, const gsl_vector* B, const gsl_vector* C, F_INC_GAMMA * figList, int debug);
double Search_F_inc(int , double , F_INC_GAMMA *);
void test(F_INC_GAMMA *f);
