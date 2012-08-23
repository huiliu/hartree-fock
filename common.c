#include "common.h"
#include <math.h>

/*
 *  放置一些常用函数
 */
void matrix_output(const gsl_matrix *m, int n, char *msg)
{
    int i, j;

    printf("%s\n",msg);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%15.6lE", gsl_matrix_get(m, i, j));
        printf("\n");
    }
}

void vector_output(const gsl_vector *v, int n, char *msg)
{
    int i;

    printf("%s\n",msg);
    for (i = 0; i < n; i++)
        printf("%15.06lE", gsl_vector_get(v, i));
    printf("\n");
}

#define F_INC_GAMMA_CYCLE    100
#define F_INC_GAMMA_delta  1.0E-10
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
        result = factorial_2(2*m -1) / pow(2*w, m + 0.5) * sqrt(M_PI_2);
    return result;
}

int factorial(int n)
{
    if (n <= 1) return 1;

    return n*factorial(n-1);
}

int factorial_2(int n)
{
    if (n % 2 == 0) {
        if (n <= 2) return 2;
    }else{
        if (n <= 1) return 1;
    }
    return n * factorial_2(n-2);
}

// check the symtery of two-electron integral
#define GetIndex(i, j, k, l, N) (N * N * N * i + N * N * j + N * k + l)
#ifndef MIN
#define MIN(a, b)   ((a) < (b) ? (a) : (b))
#endif

int ChkERISym(double ****e, int i, int j, int k, int l, int N, int *is_dup)
{
// thanks David pulq for his method
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
