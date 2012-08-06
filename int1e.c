#include <stdio.h>
#include <stdlib.h>
#include "gsl/gsl_blas.h"
#include "overlap.h"
#include "common.h"

double RecCoeff(int i, int j, int t, double *gamma, double *PA, double *PB)
{
// L. E. McMurchie, R. Davidson J. comp. pyhy. 26, 218 (1978)
// (2.20)   (2.21)
    double item1, item2, item3;
    if (i == 0 && j == 0 && t == 0) return 1.0;
    if (i > 0) {
        if (t < 0) {
            item1 = item2 = 0;
            if (t+1 <= 0)
                return 0;
            else
                return (t+1) * RecCoeff(i-1, j, t+1, gamma, PA, PB);
        }else if (t+1 > i+j) {
            //item3 = 0;
            if (t > i+j || (*PA) == 0)
                item2 = 0;
            else
                item2 = (*PA)*RecCoeff(i-1, j, t, gamma, PA, PB);

            if (t - 1 > i + j)
                return 0;
            else
                item1 = 0.5 / (*gamma) * RecCoeff(i-1, j, t-1, gamma, PA, PB);

            return item2 + item1;
        }

        if ((*PA) == 0) 
            item2 = 0;
        else
            item2 = (*PA)*RecCoeff(i-1, j, t,  gamma,PA, PB);

        item1 = 0.5/(*gamma) * RecCoeff(i-1, j, t-1, gamma, PA, PB);
        item3 = (t+1) * RecCoeff(i-1, j, t+1, gamma, PA, PB);
        return item1 + item2 + item3;
    }else if (j > 0) {
        if (t < 0) {
            item1 = item2 = 0;
            if (t+1 <= 0)
                return 0;
            else
                return (t+1) * RecCoeff(i, j-1, t+1, gamma, PA, PB);
        }else if (t+1 > i+j) {
            //item3 = 0;
            if (t > i+j || (*PB) == 0)
                item2 = 0;
            else
                item2 = (*PB)*RecCoeff(i, j-1, t, gamma, PA, PB);
            if (t - 1 > i + j)
                return 0;
            else
                item1 = 0.5/(*gamma) * RecCoeff(i, j-1, t-1, gamma, PA, PB);
            return item2 + item1;
        }

        if ((*PB) == 0) 
            item2 = 0;
        else
            item2 = (*PB)*RecCoeff(i, j-1, t, gamma, PA, PB);

        item1 = 0.5/(*gamma) * RecCoeff(i, j-1, t-1, gamma, PA, PB);
        item3 = (t+1) * RecCoeff(i, j-1, t+1, gamma, PA, PB);

        return item1 + item2 + item3;
    }
    // i == 0 && j == 0 && t < 0
    return 0;
}

double overlap_gto(const GTO* g1, const gsl_vector* A, const GTO* g2, const gsl_vector* B, int debug)
{
    double KAB, gamma;
    gsl_vector *PA, *PB;
    double result = 0;
    int l1, m1, n1, l2, m2, n2;
    double norm1, norm2, c1, c2;

    gamma = g1->alpha + g2->alpha;

    l1 = g1->l; m1 = g1->m; n1 = g1->n; norm1 = g1->norm; c1 = g1->coeff;
    l2 = g2->l; m2 = g2->m; n2 = g2->n; norm2 = g2->norm; c2 = g2->coeff;

    KAB = gauss_K(g1->alpha, A, g2->alpha, B);

    PA = gaussian_product_center(g1->alpha, A, g2->alpha, B, debug);
    PB = gsl_vector_alloc(3);
    gsl_vector_memcpy(PB, PA);
    gsl_vector_sub(PA, A);
    gsl_vector_sub(PB, B);

    result = RecCoeff(l1, l2, 0, &gamma, PA->data,   PB->data) * \
             RecCoeff(m1, m2, 0, &gamma, PA->data+1, PB->data+1) * \
             RecCoeff(n1, n2, 0, &gamma, PA->data+2, PB->data+2) * pow(M_PI/gamma, 1.5) \
             * norm1 * norm2 * c1 * c2 * KAB;
    gsl_vector_free(PA);
    gsl_vector_free(PB);
    return result;
}

double R(int n, int t, int u, int v, const gsl_vector *PX, double gamma, int debug)
{
// L. E. McMurchie, R. Davidson J. comp. pyhy. 26, 218 (1978)
// (4.6) (4.7) (4.8)
    double norm_2 = 0;
    double item1 = 0, item2 = 0;

    if (t < 0 || u < 0 || v < 0)    return 0;
    if (t == 0 && u == 0 && v == 0) {
        norm_2 =  gsl_pow_2(gsl_blas_dnrm2(PX)) * gamma;
        //if (norm_2 < 1.0E-8) norm_2 = 0;
        return gsl_pow_int(-2*gamma, n) * F_inc_gamma(n, norm_2);
    }
    if (t > 0) {
        if (t-1 == 0)
            item1 = 0;
        else
            item1 = (t-1) * R(n+1, t-2, u, v, PX, gamma, debug);
        if (PX->data[0] == 0)
            item2 = 0;
        else
            item2 = PX->data[0] * R(n+1, t-1, u, v, PX, gamma, debug);
    }
    if (u > 0) {
        if (u - 1 == 0)
            item1 = 0;
        else
            item1 = (u-1) * R(n+1, t, u-2, v, PX, gamma, debug);
        if (PX->data[1] == 0)
            item2 = 0;
        else
            item2 = PX->data[1] * R(n+1, t, u-1, v, PX, gamma, debug);
    }
    if (v > 0) {
        if (v-1 == 0)
            item1 = 0;
        else
            item1 = (v-1) * R(n+1, t, u, v-2, PX, gamma, debug);
        if (PX->data[2] == 0)
            item2 = 0;
        else
            item2 = PX->data[2] * R(n+1, t, u, v-1, PX, gamma, debug);
    }
    return item1 + item2;
}

double nuclear_elect_attraction_gto(const GTO* g1, const gsl_vector* A, \
        const GTO* g2, const gsl_vector* B, const gsl_vector* C, int debug)
{
    int i, j, k;
    gsl_vector *PC, *PA, *PB;
    double gamma, KAB, norm1, norm2, c1, c2;
    double Ex, Ey, Ez;
    int L, M, N;
    int l1, m1, n1, l2, m2, n2;
    double result = 0;

    l1 = g1->l; m1 = g1->m; n1 = g1->n; norm1 = g1->norm; c1 = g1->coeff;
    l2 = g2->l; m2 = g2->m; n2 = g2->n; norm2 = g2->norm; c2 = g2->coeff;
    
    L = l1+l2;
    M = m1+m2;
    N = n1+n2;

    gamma = g1->alpha + g2->alpha;

    KAB = gauss_K(g1->alpha, A, g2->alpha, B);

    PC = gaussian_product_center(g1->alpha, A, g2->alpha, B, debug);

    PA = gsl_vector_alloc(3);
    PB = gsl_vector_alloc(3);
    gsl_vector_memcpy(PA, PC);
    gsl_vector_memcpy(PB, PC);

    gsl_vector_sub(PA, A);
    gsl_vector_sub(PB, B);
    gsl_vector_sub(PC, C);

    for (i = 0; i <= L; i++) {
        Ex = RecCoeff(l1, l2, i, &gamma, PA->data, PB->data);
        for (j = 0; j <= M; j++) {
            Ey = RecCoeff(m1, m2, j, &gamma, PA->data+1, PB->data+1);
            for (k = 0; k <= N; k++) {
                Ez = RecCoeff(n1, n2, k, &gamma, PA->data+2, PB->data+2);
                result += Ex * Ey * Ez * R(0, i, j, k, PC, gamma,debug);
            }
        }
    }
    gsl_vector_free(PA);
    gsl_vector_free(PB);
    gsl_vector_free(PC);
    return result * norm1 * norm2 * c1 * c2 * KAB * 2 * M_PI /gamma;
}
