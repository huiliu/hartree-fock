#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>
double K_OS(double a, double b, double AB_2);
double ERI_basis_OS(const BASIS* b1, const BASIS* b2,
                    const BASIS* b3, const BASIS* b4, int debug);
double ERI_VRR_OS(int l1, int m1, int n1,
                  int l2, int m2, int n2,
                  int l3, int m3, int n3,
                  int l4, int m4, int n4,
                  double zeta, double gamma, double ro,
                  const gsl_vector *PA, const gsl_vector *PB, const gsl_vector *QC,
                  const gsl_vector *QD, const gsl_vector *WQ, const gsl_vector *WP,
                  int m, double T);
double ERI_HRR_OS(int l1, int m1, int n1,
                  int l2, int m2, int n2,
                  int l3, int m3, int n3,
                  int l4, int m4, int n4,
                  double zeta, double gamma, double ro,
                  const gsl_vector *AB, const gsl_vector *CD,
                  const gsl_vector *PA, const gsl_vector *PB, const gsl_vector *QC,
                  const gsl_vector *QD, const gsl_vector *WQ, const gsl_vector *WP,
                  double T);
