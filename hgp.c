#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "overlap.h"
#include "common.h"
#include "basis.h"
#include "hgp.h"
#include "eri_os.h"

// according contracted shell number give the number. 2*n + 1
#define MAXSHELL    11

#define GetShellBasisCount(i, L)   \
    if (i == 0) L = 1; /* S */\
    else if (i == 1) L = 3; /* P */\
    else if (i == 2) L = 6; /* D */\
    else if (i == 3) L = 10; /* F */\
    else if (i == 4) L = 15; /* G */\
    else if (i == 5) L = 21; /* H */\
    else{ \
        fprintf(stderr, "This Program cann't deal with the orbital of L >5\n"); \
        exit(EXIT_FAILURE); \
    }
        
void HGPShell(double ****ERI, const BASIS *b,
                int *ii, int *jj, int *kk, int *ll, int bc, int debug)
{
/*
 */
    int i, j, k, l;
    int L1, L2, L3, L4;
    int is_dup = 0;
    gsl_vector *AB, *CD, *B, *D;
    BASIS b1, b2, b3, b4;
    int basisCount = bc - 1;
    double *XSXS;

    AB = gsl_vector_alloc(3);
    CD = gsl_vector_alloc(3);

    B = b[*jj].xyz;
    D = b[*ll].xyz;

    gsl_vector_memcpy(AB, b[*ii].xyz);
    gsl_vector_memcpy(CD, b[*kk].xyz);

    gsl_vector_sub(AB, B);
    gsl_vector_sub(CD, D);

    // get the orbital type by angular momentum number
    GetShellBasisCount(b[*ii].L, L1);
    GetShellBasisCount(b[*jj].L, L2);
    GetShellBasisCount(b[*kk].L, L3);
    GetShellBasisCount(b[*ll].L, L4);

    MALLOC(XSXS, sizeof(double)*gsl_pow_6(MAXSHELL));

    //  use restricte condition (ab|cd)
    //if (L1 != 0 && L2 != 0 && L3 != 0 && L4 != 0) {

        for (i = 0; i < L1; i++) {
            for (j = 0; j < L2; j++) {
                for (k = 0; k < L3; k++) {
                    for (l = 0; l < L4; l++) {

                        //ChkERISym(ERI, *ii + i, *jj + j, *kk + k, *ll + l,
                        //                                        bc, &is_dup);
                        //if (is_dup)     continue;

                        b1 = b[*ii + i];
                        b2 = b[*jj + j];
                        b3 = b[*kk + k];
                        b4 = b[*ll + l];

                        if (*ll + l == 5) {
    fprintf(stdout, "%d %d %d %d\n", *ii, *jj, *kk, *ll);
        fprintf(stdout, "------------------------------------\n");
    vector_output(AB, 3, "AB:");
    vector_output(CD, 3, "CD:");
        basis_set_output(&b1, 1, "B1:");
        basis_set_output(&b2, 1, "B2:");
        basis_set_output(&b3, 1, "B3:");
        basis_set_output(&b4, 1, "B4:");
        fprintf(stdout, "------------------------------------\n");
                        }

                        ERI[*ii + i][*jj + j][*kk + k][*ll + l] = \
                           HGPBasisHRR(&b1, &b2, &b3, &b4, AB, CD, XSXS, debug);
                        
                    }
                }
            }
        }
    //}
    *ll = *ll + L4 - 1;

    if (*ll == basisCount)
        *kk = *kk + L3 - 1;
    if (*kk == basisCount)
        *jj = *jj + L2 - 1;
    if (*jj == basisCount)
        *ii = *ii + L1 - 1;

    gsl_vector_free(AB);
    gsl_vector_free(CD);
    free(XSXS);
}

double HGPBasisHRR(BASIS *b1, BASIS *b2, BASIS *b3, BASIS *b4, 
                 gsl_vector *AB, gsl_vector *CD, double *XSXS, int debug)
{
/* the input basis may be:
 *  (px, py|dx^2, dxy), ......
 * transpose it to the form of (e0|f0) by HRR(J. Chem. Phys. 89(9), 5777, 1988)
 * formula (18)
 */

// contract basis to form (e0|f0)
// (a(b+1)|cd) = ((a+1)b|cd) + AB(ab|cd)

    double item;
    if (b2->l > 0) {
        b2->l--;
        item = AB->data[0] * HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug);
        b1->l++;
        return HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug) + item;
    }
    if (b2->m > 0) {
        b2->m--;
        item = AB->data[1] * HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug);
        b1->m++;
        return HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug) + item;
    }
    if (b2->n > 0) {
        b2->n--;
        item = AB->data[1] * HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug);
        b1->n++;
        return HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug) + item;
    }

    if (b4->l > 0) {
        b4->l--;
        item = AB->data[0] * HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug);
        b3->l++;
        return HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug) + item;
    }
    if (b4->m > 0) {
        b4->m--;
        item = AB->data[1] * HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug);
        b3->m++;
        return HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug) + item;
    }
    if (b4->n > 0) {
        b4->n--;
        item = AB->data[1] * HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug);
        b3->n++;
        return HGPBasisHRR(b1, b2, b3, b4, AB, CD, XSXS, debug) + item;
    }

    // the basis integral convert to the form of (e0|f0), and continue to carry out
    // transpose from (e0|f0) to [e0|f0]
    return HGPBasis(b1, b2, b3, b4, AB, CD, XSXS, debug);
}

double HGPBasis(const BASIS* b1, const BASIS* b2,
                const BASIS* b3, const BASIS* b4,
                const gsl_vector *AB, const gsl_vector *CD,
                double *XSXS, int debug)
{
// FORM (e0|f0) FROM ∑∑∑∑[e0|f0]
    int i, j, k, l;
    int l1, m1, n1, l3, m3, n3;
    int gaussCount_1, gaussCount_2, gaussCount_3, gaussCount_4;
    int L, f;
    GTO *g1, *g2, *g3, *g4;
    double gamma1, gamma2, rho;
    gsl_vector *PA, *PB, *QC, *QD, *WQ, *WP, *PQ; 
    gsl_vector *A, *B, *C, *D;
    double norm_AB_2, norm_CD_2, T;
    double KAB, KCD;
    double *F;
    double result = 0;

    l1 = b1->l; m1 = b1->m; n1 = b1->n;
    l3 = b3->l; m3 = b3->m; n3 = b3->n;

    PB = gsl_vector_alloc(3);
    QD = gsl_vector_alloc(3);
    PQ = gsl_vector_alloc(3);

    WP = gsl_vector_alloc(3);
    WQ = gsl_vector_alloc(3);

    A = b1->xyz,
    B = b2->xyz,
    C = b3->xyz,
    D = b4->xyz,

    norm_AB_2 = gsl_pow_2(gsl_blas_dnrm2(AB));
    norm_CD_2 = gsl_pow_2(gsl_blas_dnrm2(CD));

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;
    gaussCount_3 = b3->gaussCount;
    gaussCount_4 = b4->gaussCount;

    L = l1 + m1 + n1 + l3 + m3 + n3;

    // store incomplete gamma function value
    F = calloc(sizeof(double), L + 1);

    // (e0|f0) = ∑∑∑∑[e0|f0]
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
                    rho = gamma1 * gamma2 / (gamma1 + gamma2);

                    PA = gaussian_product_center(g1->alpha, A, 
                                                 g2->alpha, B, debug);
                    QC = gaussian_product_center(g3->alpha, C,
                                                 g4->alpha, D, debug);

                    gsl_vector_memcpy(PB, PA);
                    gsl_vector_memcpy(QD, QC);
                    gsl_vector_memcpy(PQ, PA);

                    gsl_vector_sub(PQ, QC);
                    gsl_vector_memcpy(WP, PQ);
                    gsl_vector_memcpy(WQ, PQ);

                    T = rho * gsl_pow_2(gsl_blas_dnrm2(PQ));

                    for (f = 0; f <= L; f++)
                        F[f] = F_inc_gamma(f, T);

                    gsl_vector_sub(PA, A);
                    gsl_vector_sub(PB, B);
                    gsl_vector_sub(QC, C);
                    gsl_vector_sub(QD, D);

                    gsl_vector_scale(WP, -gamma2/(gamma1+gamma2));
                    gsl_vector_scale(WQ, gamma1/(gamma1+gamma2));

                    KAB = K_OS(g1->alpha, g2->alpha, norm_AB_2);
                    KCD = K_OS(g3->alpha, g4->alpha, norm_CD_2);

                    double pre = g1->norm * g1->coeff * \
                                 g2->norm * g2->coeff * \
                                 g3->norm * g3->coeff * \
                                 g4->norm * g4->coeff;

                    result += pre * KAB * KCD / sqrt(gamma1 + gamma2)* \
                                HGPHrrVRR(l1, m1, n1, l3, m3, n3,
                                          gamma1, gamma2, rho,
                                          PA, PB, QC, QD, WQ, WP, 0, F);
                    gsl_vector_free(PA);
                    gsl_vector_free(QC);
                }
            }
        }
    }
    int index = HGPIndex(l1, m1, n1, l3, m3, n3);
    XSXS[index] = result;
    free(F);
    gsl_vector_free(PB);
    gsl_vector_free(QD);
    gsl_vector_free(PQ);
    gsl_vector_free(WP);
    gsl_vector_free(WQ);

    return result;
}
int HGPIndex(l1, m1, n1, l3, m3, n3)
{
    return gsl_pow_5(MAXSHELL) * l1 + \
                gsl_pow_4(MAXSHELL)*m1 + gsl_pow_3(MAXSHELL) * n1 + \
                gsl_pow_2(MAXSHELL) * l3 + MAXSHELL * m3 + n3;
}

double HGPHrrVRR(int l1, int m1, int n1, int l3, int m3, int n3,
            double zeta, double gamma, double ro,
            const gsl_vector *PA, const gsl_vector *PB, const gsl_vector *QC,
            const gsl_vector *QD, const gsl_vector *WP, const gsl_vector *WQ,
            int m, double *T)
{
// compute the primitive integrals [e0|f0]
    double item1, item2, item3;
    double result;
    if (n3 >= 1) {

        item1 = QC->data[2] * HGPHrrVRR(l1, m1, n1, l3, m3, n3-1,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[2] * HGPHrrVRR(l1, m1, n1, l3, m3, n3-1,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n3 >= 2) {
            item2 = (n3-1)/(2*gamma) * (HGPHrrVRR(l1, m1, n1, l3, m3, n3-2,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro/gamma * HGPHrrVRR(l1, m1, n1, l3, m3, n3-2,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n1 >= 1)
            item3 = n1/(2*(zeta + gamma))*HGPHrrVRR(l1, m1, n1-1, l3, m3, n3-1,
                        zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }
    if (m3 >= 1) {

        item1 = QC->data[1] * HGPHrrVRR(l1, m1, n1, l3, m3-1, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[1] * HGPHrrVRR(l1, m1, n1, l3, m3-1, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m3 >= 2) {
            item2 = (m3-1)/(2*gamma) * (HGPHrrVRR(l1, m1, n1, l3, m3-2, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / gamma * HGPHrrVRR(l1, m1, n1, l3, m3-2, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m1 >= 1)
            item3 = m1/(2*(zeta + gamma))*HGPHrrVRR(l1, m1-1, n1, l3, m3-1, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }
    if (l3 >= 1) {

        item1 = QC->data[0] * HGPHrrVRR(l1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[0] * HGPHrrVRR(l1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l3 >= 2) {
            item2 = (l3-1) / (2*gamma) * (HGPHrrVRR(l1, m1, n1, l3-2, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * HGPHrrVRR(l1, m1, n1, l3-2, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l1 >= 1)
            item3 = l1/(2*(zeta + gamma))*HGPHrrVRR(l1-1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }

    if (n1 >= 1) {

        item1 = PA->data[2] * HGPHrrVRR(l1, m1, n1-1, l3, m3, n3,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[2] * HGPHrrVRR(l1, m1, n1-1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n1 >= 2) {
            item2 = (n1-1) / (2*zeta) * (HGPHrrVRR(l1, m1, n1-2, l3, m3, n3,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * HGPHrrVRR(l1, m1, n1-2, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n3 >= 1)
            item3 = n3/(2*(zeta + gamma))*HGPHrrVRR(l1, m1, n1-1, l3, m3, n3-1,
                        zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }
    if (m1 >= 1) {

        item1 = PA->data[1] * HGPHrrVRR(l1, m1-1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[1] * HGPHrrVRR(l1, m1-1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m1 >= 2) {
            item2 = (m1-1) / (2*zeta) * (HGPHrrVRR(l1, m1-2, n1, l3, m3, n3,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * HGPHrrVRR(l1, m1-2, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m3 >= 1)
            item3 = m3/(2*(zeta + gamma))*HGPHrrVRR(l1, m1-1, n1, l3, m3-1, n3,
                        zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }
    if (l1 >= 1) {

        item1 = PA->data[0] * HGPHrrVRR(l1-1, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[0] * HGPHrrVRR(l1-1, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l1 >= 2) {
            item2 = (l1-1) / (2*zeta) * (HGPHrrVRR(l1-2, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / zeta * HGPHrrVRR(l1-2, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l3 >= 1)
            item3 = l3/(2*(zeta + gamma))*HGPHrrVRR(l1-1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }

    return T[m];
}
