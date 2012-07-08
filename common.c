#include "common.h"

/*
 *  放置一些常用函数
 */
void matrix_output(const gsl_matrix *m, int n, char *msg)
{
    int i, j;

    printf("%s",msg);
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++)
            printf("%12.6g", gsl_matrix_get(m, i, j));
        printf("\n");
    }
}

void vector_output(const gsl_vector *v, int n, char *msg)
{
    int i;

    printf("%s",msg);
    for (i = 0; i < n; i++)
        printf("%9.06g", gsl_vector_get(v, i));
    printf("\n");
}
