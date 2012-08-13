#include <stdio.h>
#include <stdlib.h>
#include "overlap.h"
#include "int1e.h"
#include "eri.h"

double ERI_gto(const GTO* g1, const gsl_vector* A, 
                              const GTO* g2, const gsl_vector* B,
                              const GTO* g3, const gsl_vector* C,
                              const GTO* g4, const gsl_vector* D,
                              FMW *fmw, int debug)
{
    int i, j, k, l, m, n;
    int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
    int L1, M1, N1;
    int L2, M2, N2;
    gsl_vector *PA, *PB, *QC, *QD, *Q, *PQ;
    double alpha1, alpha2, alpha3, alpha4, gamma1, gamma2, finc;
    double E1x[10], E1y[10], E1z[10], E2x[10], E2y[10], E2z[10];
    double RR[20][20][20];
    double KAB, KCD;
    double result = 0;

    alpha1 = g1->alpha;
    alpha2 = g2->alpha;
    alpha3 = g3->alpha;
    alpha4 = g4->alpha;

    gamma1 = alpha1 + alpha2;
    gamma2 = alpha3 + alpha4;

    l1 = g1->l; l2 = g2->l; l3 = g3->l; l4 = g4->l;
    m1 = g1->m; m2 = g2->m; m3 = g3->m; m4 = g4->m;
    n1 = g1->n; n2 = g2->n; n3 = g3->n; n4 = g4->n;

    L1 = l1 + l2;   L2 = l3 + l4;
    M1 = m1 + m2;   M2 = m3 + m4;
    N1 = n1 + n2;   N2 = n3 + n4;


    KAB = gauss_K(alpha1, A, alpha2, B);
    KCD = gauss_K(alpha3, C, alpha4, D);

    PQ = gaussian_product_center(alpha1, A, alpha2, B, debug);
    Q = gaussian_product_center(alpha3, C, alpha4, D, debug);

    PA = gsl_vector_alloc(3); PB = gsl_vector_alloc(3);
    QC = gsl_vector_alloc(3); QD = gsl_vector_alloc(3);

    gsl_vector_memcpy(PA, PQ); gsl_vector_memcpy(PB, PQ);
    gsl_vector_memcpy(QC, Q); gsl_vector_memcpy(QD, Q);

    gsl_vector_sub(PA, A); gsl_vector_sub(PB, B);
    gsl_vector_sub(QC, C); gsl_vector_sub(QD, D);
    gsl_vector_sub(PQ, Q);

    finc = gamma1 * gamma2 / (gamma1 + gamma2);
        
    for (i = 0; i <= L1; i++)
        E1x[i] = RecCoeff(l1, l2, i, &gamma1, PA->data, PB->data);
    for (j = 0; j <= M1; j++)
        E1y[j] = RecCoeff(m1, m2, j, &gamma1, PA->data+1, PB->data+1);
    for (k = 0; k <= N1; k++)
        E1z[k] = RecCoeff(n1, n2, k, &gamma1, PA->data+2, PB->data+2);

    for (l = 0; l <= L2; l++)
        E2x[l] = RecCoeff(l3, l4, l, &gamma2, QC->data, QD->data);
    for (m = 0; m <= M2; m++)
        E2y[m] = RecCoeff(m3, m4, m, &gamma2, QC->data+1, QD->data+1);
    for (n = 0; n <= N2; n++)
        E2z[n] = RecCoeff(n3, n4, n, &gamma2, QC->data+2, QD->data+2);

    for (i = 0; i <= L1+L2; i++) {
        for (j = 0; j <= M1+M2; j++) {
            for (k = 0; k <= N1+N2; k++) {
                RR[i][j][k] = R(0, i, j, k, PQ, finc, fmw, debug); 
            }
        }
    }
    for (i = 0; i <= L1; i++) {
        for (j = 0; j <= M1; j++) {
            for (k = 0; k <= N1; k++) {
    for (l = 0; l <= L2; l++) {
        for (m = 0; m <= M2; m++) {
            for (n = 0; n <= N2; n++) {
                result += gsl_pow_int(-1, l+m+n) * E1x[i] * E1y[j] * E1z[k] * E2x[l] * E2y[m] * E2z[n] * RR[i+l][j+m][k+n];
            }
        }
    }
            }
        }
    }

    gsl_vector_free(PQ);
    gsl_vector_free(PA);
    gsl_vector_free(PB);
    gsl_vector_free(QC);
    gsl_vector_free(QD);
    gsl_vector_free(Q);

    return result * KAB * KCD *g1->norm * g2->norm * g3->norm * g4->norm * \
                    g1->coeff * g2->coeff * g3->coeff * g4->coeff \
                    * 2 * pow(M_PI, 2.5) / (gamma1 * gamma2 * sqrt(gamma1 + gamma2));
}
