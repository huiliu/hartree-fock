#include "print.h"
#include <math.h>

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

void int2e_output(double**** e, int n, char* msg)
{
    int i, j, k, l;

    printf("%s\n", msg);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++) {
                    if (fabs(e[i][j][k][l]) > 1.0E-8)
                        printf("%3d%3d%3d%3d%15.9lf\n", i, j, k, l, e[i][j][k][l]);
                }
            }
        }
    }
}

// 参数 count 表示一个基函数由count个gaussian函数组成
void basis_set_output(const BASIS* b, int count, char* msg)
{
    int i;
    printf("%s\n", msg);
    for (i = 0; i < count; i++) {
        vector_output(b[i].xyz, 3, "基组坐标:");
        gto_output(b[i].gaussian, b[i].gaussCount, "基函数");
    }
}

void gto_output(const GTO* g, int count, char* msg)
{
    int i;
    double alpha, coeff;
    int l, m, n;

    printf("%s\n", msg);
    for (i = 0; i < count; i++)
    {
        l = (g+i)->l;
        m = (g+i)->m;
        n = (g+i)->n;
        alpha = (g+i)->alpha;
        coeff = (g+i)->coeff;
        printf("%d %d %d %15.9E %16.9E\n",
                l, m, n, alpha, coeff);
    }
}

void atom_output(const ATOM_INFO** atom, int n)
{
    ATOM_INFO *a=NULL;
    int i;
    for (i = 0; i < n; i++) {
        a = (ATOM_INFO*)atom[i];
        printf("%s %d %d%12.7lf%12.7lf%12.7lf\n", a->symbol, a->n, 
                    a->basisCount, a->coordination->data[0], 
                    a->coordination->data[1], a->coordination->data[2]); 
    }
}
