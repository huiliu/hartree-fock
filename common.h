#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "basis.h"
#include <stdlib.h>

#ifndef __INTEGRAL__COMMON__
#define __INTEGRAL__COMMON__

/* -------------------------
store F_inc_gamma(m, w)
-------------------------*/
/*
typedef struct {
    int count;
    FM *FM;
}F_INC_GAMMA;

typedef struct {
    int count;
    FMW *fmw;
}FM;

typedef struct {
    double w;
    double value;
}FMW;
*/
typedef struct {
    int m;
    double w;
    double value;
}FMW;
typedef struct {
    int count;
    FMW *F;
}F_INC_GAMMA;

// 以更好的格式输出矩阵
void matrix_output(const gsl_matrix *, int , char *);
// 以更好的格式输出向量
void vector_output(const gsl_vector *, int , char *);
double F_inc_gamma(int m ,double w);
int factorial(int);
int factorial_2(int);
void* Malloc(size_t n);
void* Calloc(size_t s, size_t n);

#define MALLOC(p,n) \
    if (!(p = malloc(n))) { \
        printf("内存分配错误！malloc!\n");\
        exit(EXIT_FAILURE);\
    }

#define CALLOC(p,n,s) \
    if (!(p = calloc(n, s))) { \
        printf("内存分配错误！calloc!\n");\
        exit(EXIT_FAILURE);\
    }

#define REALLOC(p,n) \
    if (!(p = realloc(p, n))) { \
        printf("内存分配错误！ Realloc!\n");\
        exit(EXIT_FAILURE);\
    }

#endif
