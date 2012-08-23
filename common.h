#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "basis.h"
#include <stdlib.h>

#ifndef __INTEGRAL__COMMON__
#define __INTEGRAL__COMMON__
// 以更好的格式输出矩阵
void matrix_output(const gsl_matrix *, int , char *);
// 以更好的格式输出向量
void vector_output(const gsl_vector *, int , char *);
double F_inc_gamma(int m ,double w);
int factorial(int);
int factorial_2(int);
int ChkERISym(double ****e, int i, int j, int k, int l, int N, int *is_dup);
void* Malloc(size_t n);
void* Calloc(size_t s, size_t n);

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
