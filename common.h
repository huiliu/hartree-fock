#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "basis.h"
#include "redblack.h"
#include <stdlib.h>

#ifndef __INTEGRAL__COMMON__
#define __INTEGRAL__COMMON__

// use rbtree store and searc the value of incomplete gamma function.
#define ERI_INT_USE_REDBLACK_TREE

typedef struct {
    unsigned int count;
    // real count in fmwv
    struct rbtree *RBList;
}FMW;

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

#define REALLOC(p,n) \
    if (!(p = realloc(p, n))) { \
        printf("内存分配错误！\n");\
        exit(EXIT_FAILURE);\
    }

// 以更好的格式输出矩阵
void matrix_output(const gsl_matrix *, int , const char *, char *);
// 以更好的格式输出向量
void vector_output(const gsl_vector *, int , const char *, char *);
double F_inc_gamma(int m ,double w, struct rbtree *rb);
int factorial(int);
int factorial_2(int);
void* Malloc(size_t n);
void* Calloc(size_t s, size_t n);
int compare(const void *pa, const void *pb, const void *config);
unsigned int cutoff(double w);
unsigned int cutoff_small(double w);
char *replace(char *src, char *a, char *b);
#endif
