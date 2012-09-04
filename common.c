#include "common.h"
#include <math.h>

/*
 *  放置一些常用函数
 */
#define M_SQRT_PI_2             1.253314137315500121
#define F_INC_GAMMA_CYCLE       100
#define F_INC_GAMMA_delta       1.0E-10
double F_inc_gamma(int m ,double w)
{
    double result = 0;
    double tmp = 0;
    int i;
    
    if (w < 17) {
        result = tmp = 1.0 / factorial_2(2*m + 1);
        for (i = 1; i < F_INC_GAMMA_CYCLE; i++) {
            tmp *= ((2*w) / (2*m + 2*i + 1));
            if ((tmp - F_INC_GAMMA_delta) < 0)
                break;
            result += tmp;
        }
        return result * factorial_2(2 * m -1) * exp(-w);;
    }else
        result = factorial_2(2*m -1) / pow(2*w, m + 0.5) * M_SQRT_PI_2;
    return result;
}

inline int factorial(int n)
{
    int i, result = 1;
    if (n <= 1) return 1;
    for (i = 2; i <= n; i++)
        result *= i;
    return result;
}

double fact2[] = {1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025, 10321920, 34459425, 185794560, 654729075, 3715891200, 13749310575, 81749606400, 316234143225, 1961990553600};

inline int factorial_2(int n)
{
    if (n <= 0) return 1;
    return fact2[n];
}

// check the symtery of two-electron integral
#define GetIndex(i, j, k, l, N) (N * N * N * i + N * N * j + N * k + l)
#ifndef MIN
#define MIN(a, b)   ((a) < (b) ? (a) : (b))
#endif

inline int ChkERISym(double ****e, int i, int j, int k, int l, int N, int *is_dup)
{
// thanks David pulq for his idea.
// https://plus.google.com/106075773891428215861/posts
// https://gist.github.com/3427265
    long pos;
    long min;

    if (e[i][j][k][l] != 0){
        *is_dup = 1;
        return 0;
    }

    pos = min = GetIndex(i, j, k, l, N);
    
    min = MIN(min, GetIndex(i, j, l, k, N));
    min = MIN(min, GetIndex(j, i, l, k, N));
    min = MIN(min, GetIndex(j, i, k, l, N));
    min = MIN(min, GetIndex(k, l, i, j, N));
    min = MIN(min, GetIndex(k, l, j, i, N));
    min = MIN(min, GetIndex(l, k, i, j, N));
    min = MIN(min, GetIndex(l, k, j, i, N));
    
    *is_dup = pos != min;
    
    return 0;
}

void* Malloc(size_t n)
{
    void * p;
    if (!(p = malloc(n))) {
        printf("Memory Assigement Wrong!\n");
        exit(EXIT_FAILURE);
    }else
        return p;
}

void* Calloc(size_t s, size_t n)
{
    void * p;
    if (!(p = calloc(s, n))) {
        printf("Memory Assigement Wrong!\n");
        exit(EXIT_FAILURE);
    }else
        return p;
}

inline gsl_vector* gaussian_product_center(const double a, const gsl_vector *A, 
                            const double b, const gsl_vector *B, int flags)
{
// Gaussian函数乘积定理计算双中心
    int i;
    double gamma = a + b;
    gsl_vector *center = gsl_vector_alloc(3);
    
    for (i = 0; i < 3; i++)
        center->data[i] = (a * A->data[i] + b * B->data[i]) / gamma;

    return center;
}
