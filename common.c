#include "common.h"
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

/*
 *  放置一些常用函数
 */
void matrix_output(const gsl_matrix *m, int n, const char *c, char *suffix)
{
    int i, j, len;
    FILE *f;
    char randchar[3], fName[50];
    time_t second = time(NULL) % 17;

    if (suffix == NULL)
        f = stdout;
    else{
        len = strlen(c);
        strcpy(fName, c);
        replace(fName, "inp", suffix);
        sprintf(randchar, "%ld", second);
        strcat(fName, randchar);

        f = fopen(fName, "w");
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            fprintf(f, "%15.6lE", gsl_matrix_get(m, i, j));
        fprintf(f, "\n");
    }
    //fclose(f);
}

void vector_output(const gsl_vector *v, int n, const char *c, char *suffix)
{
    int i, j, len;
    FILE *f;
    char randchar[3], fName[50];
    time_t second = time(NULL) % 17;

    if (suffix == NULL)
        f = stdout;
    else {
        len = strlen(c);
        strcpy(fName, c);
        replace(fName, "inp", suffix);
        sprintf(randchar, "%ld", second);
        strcat(fName, randchar);

        f = fopen(fName, "w");
    }

    for (i = 0; i < n; i++)
        fprintf(f, "%15.06lE", gsl_vector_get(v, i));
    //fclose(f);
}

#define F_INC_GAMMA_CYCLE    100
#define F_INC_GAMMA_delta  1.0E-10
double F_inc_gamma(int m ,double w, struct rbtree *rb)
{

    double result = 0;
#ifdef ERI_INT_USE_REDBLACK_TREE
    ITEM *f_result;
    MALLOC(f_result, sizeof(ITEM));

    f_result->w = w;
    f_result->value =  -999;
    // search in the redblack tree
    f_result = (ITEM *)rbsearch(f_result, rb);

    if (f_result->value >= 0)
        return f_result->value;
#endif
    double tmp = 0;
    register int i;
    if (w <= 17) {
        if (m == 0 && w == 0)  result = 1.0;

        result = tmp = 1.0 / factorial_2(2*m + 1);
        for (i = 1; i < F_INC_GAMMA_CYCLE; i++) {
            tmp *= ((2*w) / (2*m + 2*i + 1));
            if ((tmp - F_INC_GAMMA_delta) < 0)
                break;
            result += tmp;
        }
        result *= factorial_2(2 * m -1) * exp(-w);;
    }else
        result = factorial_2(2*m -1) / pow(2*w, m + 0.5) * sqrt(M_PI_2);

#ifdef ERI_INT_USE_REDBLACK_TREE
    f_result->value = result;
#endif

    return result;
}

int factorial(int n)
{
    if (n <= 1) return 1;
    register int i = n;
    register int result = 1;

    for (; i > 0; i--)
        result *= i;
    return result;
}

int factorial_2(int n)
{
    register int i = n;
    register int result = 1;

    if (i % 2 == 0) {
        if (i <= 2) return 2;
    }else{
        if (i <= 1) return 1;
    }

    for (; i > 0; i -= 2) {
        result *= i;
    }
    return result;
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

int compare(const void *pa, const void *pb, const void *config)
{
#define epsilon 1.0E-8
    ITEM *a, *b;
    a = (ITEM *)pa; b = (ITEM *)pb;
    return gsl_fcmp(a->w, b->w, epsilon);
	if(a->w < b->w) return -1;
	if(a->w > b->w) return 1;
    return 0;
}

char *replace(char *src, char *a, char *b)
{
    assert(src); assert(a); assert(b);

    char *tmp;
    int i = strlen(a);

    tmp = strstr(src, a);
    strncpy(tmp, b, i);
    return src;
}
