#include "common.h"
#include "overlap.h"
#include "basis.h"


#ifndef __INTEGRAL__INT2E__
#define __INTEGRAL__INT2E__

double theta(int l, int l1, int l2, double PA, double PB, int r, double gamma);
double B(int l, int l1, int l2, double PA, double PB, int r, double gamma1,
         int ll,int l3, int l4, double QC, double QD, int rr,double gamma2,
         int i, double px);
double omega(double alpha1, double alpha2, double alpha3, double alpha4, 
                                                        double AB, double CD);
void Bxyz(int l1, int l2, double PA, double PB, double gamma1, 
          int l3, int l4, double QC, double QD, double gamma2,
          double px, double* result);
double int2e_gto(const GTO* g1, const gsl_vector* A, 
                              const GTO* g2, const gsl_vector* B,
                              const GTO* g3, const gsl_vector* C,
                              const GTO* g4, const gsl_vector* D, int debug);
double int2e_basis(const BASIS* b1, const BASIS* b2,
                                const BASIS* b3, const BASIS* b4);
double* int2e_matrix(INPUT_INFO* b);
void int2e_output(double* e, int n, char* msg);
#endif
