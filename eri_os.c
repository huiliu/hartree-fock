#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "common.h"
#include "basis.h"
#include "eri_os.h"
#include "print.h"

double ERI_basis_OS(const BASIS* b1, const BASIS* b2,
                    const BASIS* b3, const BASIS* b4, COORD **PP, int debug)
{ 
    int i, j, k, l, gaussCount_1, gaussCount_2, gaussCount_3, gaussCount_4;
    int L, f;
    GTO *g1, *g2, *g3, *g4;
    double gamma1, gamma2, ro;
    gsl_vector *P, *Q, *PA, *PB, *QC, *QD, *WQ, *WP, *PQ; 
    gsl_vector *A, *B, *C, *D, *AB, *CD;
    double norm_PQ_2, norm_AB_2, norm_CD_2, T;
    double KAB, KCD, pre;
    double *F;
    double result = 0;

    P = gsl_vector_alloc(3);
    Q = gsl_vector_alloc(3);
    PB = gsl_vector_alloc(3);
    QD = gsl_vector_alloc(3);
    PQ = gsl_vector_alloc(3);
    AB = gsl_vector_alloc(3);
    CD = gsl_vector_alloc(3);

    WP = gsl_vector_alloc(3);
    WQ = gsl_vector_alloc(3);

    A = b1->xyz,
    B = b2->xyz,
    C = b3->xyz,
    D = b4->xyz,

/*
    basis_set_output(b1, 1, "A:");
    basis_set_output(b2, 1, "B:");
    basis_set_output(b3, 1, "C:");
    basis_set_output(b4, 1, "D:");
*/

    gsl_vector_memcpy(AB, A);
    gsl_vector_sub(AB, B);
    gsl_vector_memcpy(CD, C);
    gsl_vector_sub(CD, D);

    norm_AB_2 = gsl_pow_2(gsl_blas_dnrm2(AB));
    norm_CD_2 = gsl_pow_2(gsl_blas_dnrm2(CD));

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;
    gaussCount_3 = b3->gaussCount;
    gaussCount_4 = b4->gaussCount;

    L = b1->l + b1->m + b1->n + b2->l + b2->m + b2->n + b3->l + b3->m + b3->n + b4->l + b4->m + b4->n;

    F = calloc(sizeof(double), L + 1);

    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            for (k = 0; k < gaussCount_3; k++) {
                for (l = 0; l < gaussCount_4; l++) {

                    g1 = &b1->gaussian[i];
                    g2 = &b2->gaussian[j];
                    g3 = &b3->gaussian[k];
                    g4 = &b4->gaussian[l];

                    gamma1 = g1->alpha + g2->alpha;
                    gamma2 = g3->alpha + g4->alpha;
                    ro = gamma1 * gamma2 / (gamma1 + gamma2);

                    PA = gaussian_product_center(g1->alpha, A, g2->alpha, B, debug);
                    QC = gaussian_product_center(g3->alpha, C, g4->alpha, D, debug);

                    gsl_vector_memcpy(P, PA);
                    gsl_vector_memcpy(Q, QC);
                    
                    gsl_vector_memcpy(PB, PA);
                    gsl_vector_memcpy(QD, QC);
                    gsl_vector_memcpy(PQ, PA);

                    gsl_vector_sub(PQ, QC);
                    gsl_vector_memcpy(WP, PQ);
                    gsl_vector_memcpy(WQ, PQ);

                    norm_PQ_2 = gsl_pow_2(gsl_blas_dnrm2(PQ));
                    T = ro * norm_PQ_2;

                    for (f = 0; f <= L; f++)
                        F[f] = F_inc_gamma(f, T);

                    //fprintf(stdout, "%16.9E%16.9E%16.9E%16.9E%16.9E\n", ro, PQ->data[0], PQ->data[1], PQ->data[2], T);

                    /*
                    F[L] = F_inc_gamma(L, T);
                    for (f = L-1; f >= 0; f--)
                        F[f] = (2 * T * F[f+1] + exp(-T)) / (2*f + 1);
                    */

                    gsl_vector_sub(PA, A);
                    gsl_vector_sub(PB, B);
                    gsl_vector_sub(QC, C);
                    gsl_vector_sub(QD, D);

                    gsl_vector_scale(WP, -gamma2/(gamma1+gamma2));
                    gsl_vector_scale(WQ, gamma1/(gamma1+gamma2));

                    KAB = K_OS(g1->alpha, g2->alpha, norm_AB_2);
                    KCD = K_OS(g3->alpha, g4->alpha, norm_CD_2);

                    pre = g1->coeff * g2->coeff * g3->coeff * g4->coeff;

/*
                    fprintf(stdout, "----- ----- ----- -----\n \
                                     %d %d %d\n \
                                     %d %d %d\n \
                                     %d %d %d\n \
                                     %d %d %d\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E\n", 
                                    g1->l, g1->m ,g1->n,
                                    g2->l, g2->m ,g2->n,
                                    g3->l, g3->m ,g3->n,
                                    g4->l, g4->m ,g4->n,
                                    gamma1, gamma2, ro,
                                    pre, KAB, KCD,
                                    P->data[0], P->data[1], P->data[2],
                                    Q->data[0], Q->data[1], Q->data[2],
                                    PA->data[0], PA->data[1], PA->data[2],
                                    PB->data[0], PB->data[1], PB->data[2],
                                    QC->data[0], QC->data[1], QC->data[2],
                                    QD->data[0], QD->data[1], QD->data[2],
                                    WQ->data[0], WQ->data[1], WQ->data[2],
                                    WP->data[0], WP->data[1], WP->data[2],
                                    T);
*/                                    
/*
fprintf(stdout, "===== ===== =====\n%15.9E %15.9E %15.9E %15.9E\n",
                                    g1->alpha, g2->alpha, g3->alpha, g4->alpha);
fprintf(stdout, "----- ----- -----\n%15.9E %15.9E\n", gamma1, gamma2);
*/


                    result += pre * KAB * KCD * ERI_VRR_OS(g1->l, g1->m ,g1->n,
                                                       g2->l, g2->m ,g2->n,
                                                       g3->l, g3->m ,g3->n,
                                                       g4->l, g4->m ,g4->n,
                                                       gamma1, gamma2, ro,
                                                       /*AB, CD,*/
                                                       PA, PB, QC, QD, WQ, WP,
                                                       0, F)/sqrt(gamma1 + gamma2);
                    gsl_vector_free(PA);
                    gsl_vector_free(QC);
                }
            }
        }
    }
    free(F);
    gsl_vector_free(P);
    gsl_vector_free(Q);
    gsl_vector_free(PB);
    gsl_vector_free(QD);
    gsl_vector_free(PQ);
    gsl_vector_free(AB);
    gsl_vector_free(CD);
    gsl_vector_free(WP);
    gsl_vector_free(WQ);
    if (debug == 999)
        fprintf(stdout, "result: %15.8lf\n", result);
    return result;
}

double ERI_VRR_OS(int l1, int m1, int n1,
                  int l2, int m2, int n2,
                  int l3, int m3, int n3,
                  int l4, int m4, int n4,
                  double zeta, double gamma, double ro,
                  const gsl_vector *PA, const gsl_vector *PB, const gsl_vector *QC,
                  const gsl_vector *QD, const gsl_vector *WQ, const gsl_vector *WP,
                  int m, double *T)
{
    double item1, item2, item3, item4, item5, result;
    // -----------------------------------the fourth GTO-------------------------------
    if (n4 >= 1) {

        item1 = QD->data[2] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[2] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n4 >= 2) {
            item2 = (n4-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-2, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-2, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n3 >=1)
            item3 = n3 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (n2 >= 1)
            item4 = n2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (n1 >= 1)
            item5 = n1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    if (m4 >= 1) {

        item1 = QD->data[1] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[1] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m4 >= 2) {
            item2 = (m4-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-2, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-2, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m3 >=1)
            item3 = m3 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (m2 >= 1)
            item4 = m2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (m1 >= 1)
            item5 = m1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    if (l4 >= 1) {

        item1 = QD->data[0] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[0] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l4 >= 2) {
            item2 = (l4-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-2, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-2, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l3 >=1)
            item3 = l3 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (l2 >= 1)
            item4 = l2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (l1 >= 1)
            item5 = l1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }
    // -----------------------------------the third GTO-------------------------------
    if (n3 >= 1) {

        item1 = QC->data[2] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[2] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n3 >= 2) {
            item2 = (n3-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-2, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-2, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n4 >=1)
            item3 = n4 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (n1 >= 1)
            item4 = n1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (n2 >= 1)
            item5 = n2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    if (m3 >= 1) {

        item1 = QC->data[1] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[1] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m3 >= 2) {
            item2 = (m3-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-2, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-2, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m4 >=1)
            item3 = m4 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (m1 >= 1)
            item5 = m1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (m2 >= 1)
            item4 = m2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    if (l3 >= 1) {

        item1 = QC->data[0] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[0] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l3 >= 2) {
            item2 = (l3-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-2, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-2, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l4 >=1)
            item3 = l4 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (l1 >= 1)
            item5 = l1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (l2 >= 1)
            item4 = l2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }
    // ---------------------------------the second basis-------------------------------------
    if (n2 >= 1) {

        item1 = PB->data[2] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[2] * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n2 >= 2) {
            item2 = (n2-1) / (2*zeta) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2-2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n1 >=1)
            item3 = n1 / (2*zeta) * (ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (n4 >= 1)
            item4 = n4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP,  m+1, T);
        else
            item4 = 0;

        if (n3 >= 1)
            item5 = n1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP,  m+1, T);
        else
            item5 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    if (m2 >= 1) {

        item1 = PB->data[1] * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[1] * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m2 >= 2) {
            item2 = (m2-1) / (2*zeta) * (ERI_VRR_OS(l1, m1, n1, l2, m2-2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1, l2, m2-2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m1 >=1)
            item3 = m1 / (2*zeta) * (ERI_VRR_OS(l1, m1-1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1-1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (m4 >= 1)
            item5 = m4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (m3 >= 1)
            item4 = m3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    if (l2 >= 1) {

        item1 = PB->data[0] * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[0] * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l2 >= 2) {
            item2 = (l2-1) / (2*zeta) * (ERI_VRR_OS(l1, m1, n1, l2-2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1, l2-2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l1 >=1)
            item3 = l1 / (2*zeta) * (ERI_VRR_OS(l1-1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1-1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (l4 >= 1)
            item5 = l4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (l3 >= 1)
            item4 = l3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    // the first basis
    if (n1 >= 1) {

        item1 = PA->data[2] * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[2] * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n1 >= 2) {
            item2 = (n1-1) / (2*zeta) * (ERI_VRR_OS(l1, m1, n1-2, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1-2, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n2 >=1)
            item3 = n2 / (2*zeta) * (ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (n3 >= 1)
            item4 = n3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4,  zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (n4 >= 1)
            item5 = n4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP,  m+1, T);
        else
            item5 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    if (m1 >= 1) {

        item1 = PA->data[1] * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[1] * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m1 >= 2) {
            item2 = (m1-1) / (2*zeta) * (ERI_VRR_OS(l1, m1-2, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1-2, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m2 >=1)
            item3 = m2 / (2*zeta) * (ERI_VRR_OS(l1, m1-1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1-1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (m3 >= 1)
            item5 = m3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (m4 >= 1)
            item4 = m4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

    if (l1 >= 1) {

        item1 = PA->data[0] * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[0] * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l1 >= 2) {
            item2 = (l1-1) / (2*zeta) * (ERI_VRR_OS(l1-2, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1-2, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l2 >=1)
            item3 = l2 / (2*zeta) * (ERI_VRR_OS(l1-1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1-1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (l3 >= 1)
            item5 = l3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (l4 >= 1)
            item4 = l4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        result = item1 + item2 + item3 + item4 + item5;
        return result;
    }

        //fprintf(stdout, "F%8d%20.8lf%20.8lf\n", m, T, F_inc_gamma(m, T));
    return T[m];
}
