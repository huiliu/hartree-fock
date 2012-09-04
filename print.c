#include "print.h"

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

#ifdef __INTEGRAL__INT2E__ONE__
void int2e_output(double* e, int n, char* msg)
#else
void int2e_output(double**** e, int n, char* msg)
#endif
{
    int i, j, k, l;
#ifdef __INTEGRAL__INT2E__ONE__
    int I;
#endif

    printf("%s\n", msg);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++) {
#ifdef __INTEGRAL__INT2E__ONE__
                    I = i * gsl_pow_3(n) + j * gsl_pow_2(n) + \
                                                            k * n + l;
                    printf("%15.9lf", e[I]);
#else
                    if (fabs(e[i][j][k][l]) > 1.0E-8)
                        printf("%3d%3d%3d%3d%15.9lf\n", i, j, k, l, e[i][j][k][l]);
#endif
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
    double alpha, coeff, norm;
    int l, m, n;

    printf("%s\n", msg);
    for (i = 0; i < count; i++)
    {
        l = (g+i)->l;
        m = (g+i)->m;
        n = (g+i)->n;
        alpha = (g+i)->alpha;
        coeff = (g+i)->coeff;
        norm = (g+i)->norm;
        printf("%d %d %d %12.8lf %12.8lf %12.8lf\n",
                l, m, n, alpha, coeff, norm);
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
        // for (j = 0; j < a->basisCount; j++)
        //     basis_set_output(a->basis + j, 3, "Basis Function:");
    }
}
