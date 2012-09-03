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

double HGPShell(const BASIS *b,
                int ii, int jj, int kk, int ll, int debug)
{
/*
 * select all basis function in common shell. for example px, py, pz;
 *
 */
    gsl_vector *AB, *CD;
    int L1, L2, L3, L4;
    int i, j, k, l;
    int r, s, u, v;
    BASIS b1, b2, b3, b4;
    double *XSXS = NULL;
    double result;

    AB = gsl_vector_alloc(3);
    CD = gsl_vector_alloc(3);

    i = j = k = l = 0;
    r = ii + i;
    s = jj + j;
    u = kk + k;
    v = ll + l;

    b1 = b[r];
    b2 = b[s];
    b3 = b[u];
    b4 = b[v];

    L1 = b1.L;
    L2 = b2.L;
    L3 = b3.L;
    L4 = b4.L;

    gsl_vector_memcpy(AB, b1.xyz);
    gsl_vector_memcpy(CD, b3.xyz);

    gsl_vector_sub(AB, b2.xyz);
    gsl_vector_sub(CD, b4.xyz);

    result = HGPBasisHRR(b1, b2, b3, b4,
                             b1.l, b1.m, b1.n,
                             b2.l, b2.m, b2.n,
                             b3.l, b3.m, b3.n,
                             b4.l, b4.m, b4.n,
                             AB, CD, XSXS, debug);
                        
    gsl_vector_free(AB);
    gsl_vector_free(CD);

    return result;
    //free(XSXS);
}

double HGPBasisHRR(BASIS b1, BASIS b2, BASIS b3, BASIS b4, 
                    int l1, int m1, int n1,
                    int l2, int m2, int n2,
                    int l3, int m3, int n3,
                    int l4, int m4, int n4,
                 gsl_vector *AB, gsl_vector *CD, double *XSXS, int debug)
{
/* the input basis may be:
 * (px, py|dx^2, dxy), ......
 * transpose it to the form of (e0|f0) by HRR(J. Chem. Phys. 89(9), 5777, 1988)
 * formula (18)
 */

// contract basis to form (e0|f0)

// (a(b+1)|cd) = ((a+1)b|cd) + AB(ab|cd)
    if (l2 > 0) {
        return          HGPBasisHRR(b1, b2, b3, b4, l1+1,m1,n1,l2-1,m2,n2,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        //+ AB->data[0] * HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2-1,m2,n2,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
    }
    if (m2 > 0) {
        return          HGPBasisHRR(b1, b2, b3, b4, l1,m1+1,n1,l2,m2-1,n2,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        //+ AB->data[1] * HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2,m2-1,n2,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
    }
    if (n2 > 0) {
        return          HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1+1,l2,m2,n2-1,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        //+ AB->data[2] * HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2,m2,n2-1,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
    }

// (ab|c(d+1)) = (ab|(c+1)d) + CD(ab|cd)
    if (l4 > 0) {
        return          HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2,m2,n2,l3+1,m3,n3,l4-1,m4,n4,AB, CD, XSXS, debug);
        //+ CD->data[0] * HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2,m2,n2,l3,m3,n3,l4-1,m4,n4,AB, CD, XSXS, debug);
    }
    if (m4 > 0) {
        return          HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2,m2,n2,l3,m3+1,n3,l4,m4-1,n4,AB, CD, XSXS, debug);
        //+ CD->data[1] * HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2,m2,n2,l3,m3,n3,l4,m4-1,n4,AB, CD, XSXS, debug);
    }
    if (n4 > 0) {
        return          HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2,m2,n2,l3,m3,n3+1,l4,m4,n4-1,AB, CD, XSXS, debug);
        //+ CD->data[2] * HGPBasisHRR(b1, b2, b3, b4, l1,m1,n1,l2,m2,n2,l3,m3,n3,l4,m4,n4-1,AB, CD, XSXS, debug);
    }

    // the basis integral has been converted to the form of (e0|f0), 
    // and continue to carry out transpose (e0|f0) to [e0|f0]
    // FORM (e0|f0) FROM ∑∑∑∑[e0|f0]
    return HGPBasis(&b1, &b2, &b3, &b4, l1, m1, n1, l3, m3, n3, AB, CD, XSXS, debug);
}

void HGPBra(BASIS *b1, BASIS *b2, BASIS *bra)
{
    ;
}

double HGPBasis(const BASIS* b1, const BASIS* b2,
                const BASIS* b3, const BASIS* b4,
                int l1, int m1, int n1, int l3, int m3, int n3,
                const gsl_vector *AB, const gsl_vector *CD,
                double *XSXS, int debug)
{
    // [e0|f0] ---> (e0|f0)
    int i, j, k, l;
    //int l1, m1, n1, l3, m3, n3;
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

    //l1 = b1->l; m1 = b1->m; n1 = b1->n;
    //l3 = b3->l; m3 = b3->m; n3 = b3->n;


    PB = gsl_vector_alloc(3);
    QD = gsl_vector_alloc(3);
    PQ = gsl_vector_alloc(3);

    WP = gsl_vector_alloc(3);
    WQ = gsl_vector_alloc(3);

    A = b1->xyz;
    B = b2->xyz;
    C = b3->xyz;
    D = b4->xyz;

    if (debug == 5) {
        vector_output(A, 3, "A:");
        vector_output(B, 3, "B:");
        vector_output(C, 3, "C:");
        vector_output(D, 3, "D:");
        vector_output(AB, 3, "AB:");
        vector_output(CD, 3, "CD:");
    }

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

        vector_output(PA, 3, "PA:");
        vector_output(PB, 3, "PB:");
        vector_output(QC, 3, "QC:");
        vector_output(QD, 3, "QD:");
        vector_output(WP, 3, "WP:");
        vector_output(WQ, 3, "WQ:");

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
    //int index = HGPIndex(l1, m1, n1, l3, m3, n3);
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
// VRR
// compute the primitive integrals [e0|f0]
    double item1 = 0, item2, item3;
    double result;
    if (n3 >= 1) {
        /*
        item1 = QC->data[2] * HGPHrrVRR(l1, m1, n1, l3, m3, n3-1,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[2] * HGPHrrVRR(l1, m1, n1, l3, m3, n3-1,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        */

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
        /*
        item1 = QC->data[1] * HGPHrrVRR(l1, m1, n1, l3, m3-1, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[1] * HGPHrrVRR(l1, m1, n1, l3, m3-1, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        */
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
        /*
        item1 = QC->data[0] * HGPHrrVRR(l1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->data[0] * HGPHrrVRR(l1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        */

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
        /*
        item1 = PA->data[2] * HGPHrrVRR(l1, m1, n1-1, l3, m3, n3,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[2] * HGPHrrVRR(l1, m1, n1-1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        */
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
        /*
        item1 = PA->data[1] * HGPHrrVRR(l1, m1-1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[1] * HGPHrrVRR(l1, m1-1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        */
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
        /*
        item1 = PA->data[0] * HGPHrrVRR(l1-1, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->data[0] * HGPHrrVRR(l1-1, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        */
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
