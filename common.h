#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "basis.h"
#include <stdlib.h>

#ifndef __INTEGRAL__COMMON__
#define __INTEGRAL__COMMON__
double F_inc_gamma(int m ,double w);
inline int factorial(int);
inline int factorial_2(int);
// check the symtery of two-electron integral
inline int ChkERISym(double ****e, int i, int j, int k, int l, int N, int *is_dup);
void* Malloc(size_t n);
void* Calloc(size_t s, size_t n);
inline gsl_vector* gaussian_product_center(const double a, const gsl_vector *A, const double b, const gsl_vector *B, int);

#define MALLOC(p,n) \
    if (!(p = malloc(n))) { \
        printf("内存分配错误！\n");\
        exit(EXIT_FAILURE);\
    }

#define CALLOC(p,n,s) \
    if (!(p = calloc(n, s))) { \
        printf("内存分配错误！\n");\
        exit(EXIT_FAILURE);\
    }

#endif
