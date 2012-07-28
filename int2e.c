#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "overlap.h"
#include "basis.h"
#include "int2e.h"
#include <gsl/gsl_blas.h>

double theta(int l, int l1, int l2, double PA, double PB, int r, double gamma)
{
// David B. Cook; Handbook of Computational Quantum Chemistry. P249
    int debug = 0;
    return fi_l_ll_pax_pbx(l, l1, l2, PA, PB, debug) * factorial(l) * \
            pow(gamma, r-l) / factorial(r) / factorial(l - 2*r);
}

double B(int l, int l1, int l2, double PA, double PB, int r, double gamma1,
         int ll,int l3, int l4, double QC, double QD, int rr,double gamma2,
         int i, double px)
{
    double result = 0;
    double delta = 0.25 * (1/gamma1 + 1/gamma2);

    result = pow(-1, ll+i) * theta(l, l1, l2, PA, PB, r, gamma1) * \
                            theta(ll, l3, l4, QC, QD, rr, gamma2);
    result *= pow(2*delta, 2*(r+rr)) * factorial(l + ll -2*(r + rr)) * \
                pow(delta, i) * pow(px, l + ll - 2*(r + rr + i));
    result = result / pow(4*delta, l + ll) * factorial(i) * factorial(l + ll - \
                                                                2*(r + rr + i));

    return result;
}

double omega(double alpha1, double alpha2, double alpha3, double alpha4, 
                                                        double AB, double CD)
{
    double gamma1 = alpha1 + alpha2;
    double gamma2 = alpha3 + alpha4;

    return 2*pow(M_PI, 2)/(gamma1*gamma2) * sqrt(M_PI / (gamma1+gamma2)) * exp(
                                        -alpha1*alpha2 * pow(AB, 2) / gamma1 
                                        -alpha3*alpha4 * pow(CD, 2) / gamma2);
}

void Bxyz(int l1, int l2, double PA, double PB, double gamma1, 
          int l3, int l4, double QC, double QD, double gamma2,
          double px, double* result)
{
    int l, r, i, ll, rr, I;


    for (l = 0; l <= l1+l2; l++) {
        for (r = 0; r <= l/2; r++) {
            for (i = 0; i <= (l - 2*r)/2; i++) {
                for (ll = 0; ll <= l3 + l4; ll++) {
                    for (rr = 0; rr <= ll/2; rr++) {
                        I = l + ll - 2*r + 2*rr - i;
                        result[I] += B(l, l1, l2, PA, PB, r, gamma1, 
                                      ll, l3, l4, QC, QD, rr, gamma2, i, px);
                    }
                }
            }
        }
    }
}

double int2e_gto(const GTO* g1, const gsl_vector* A, 
                              const GTO* g2, const gsl_vector* B,
                              const GTO* g3, const gsl_vector* C,
                              const GTO* g4, const gsl_vector* D, int debug)
{
    int i, j, k;
    int L, M, N;
    gsl_vector *PA, *PB, *QC, *QD, *P, *Q, *PQ, *AB, *CD;
    double Bx[10], By[10], Bz[10];
    double alpha1, alpha2, alpha3, alpha4, gamma1, gamma2;
    double norm_pq_2, norm_ab, norm_cd, finc;
    double result = 0;

    alpha1 = g1->alpha;
    alpha2 = g2->alpha;
    alpha3 = g3->alpha;
    alpha4 = g4->alpha;

    gamma1 = alpha1 + alpha2;
    gamma2 = alpha3 + alpha4;

    L = g1->l + g2->l + g3->l + g4->l;
    M = g1->m + g2->m + g3->m + g4->m;
    N = g1->n + g2->n + g3->n + g4->n;

    P = gaussian_product_center(g1->alpha, A, g2->alpha, B, debug);
    Q = gaussian_product_center(g3->alpha, C, g4->alpha, D, debug);

    PA = gsl_vector_alloc(3); PB = gsl_vector_alloc(3);
    QC = gsl_vector_alloc(3); QD = gsl_vector_alloc(3);
    AB = gsl_vector_alloc(3); CD = gsl_vector_alloc(3);

    gsl_vector_memcpy(PA, P); gsl_vector_memcpy(PB, P);
    gsl_vector_memcpy(QC, Q); gsl_vector_memcpy(QD, Q);

    PQ = gsl_vector_alloc(3);
    gsl_vector_memcpy(PQ, P);
    gsl_vector_sub(PQ, Q);
    norm_pq_2 = pow(gsl_blas_dnrm2(PQ), 2);

    gsl_vector_sub(PA, A); gsl_vector_sub(PB, B);
    gsl_vector_sub(QC, C); gsl_vector_sub(QD, D);

    gsl_vector_memcpy(AB, A); gsl_vector_memcpy(CD, C);
    gsl_vector_sub(AB, B); gsl_vector_sub(CD, D);

    norm_ab = gsl_blas_dnrm2(AB);
    norm_cd = gsl_blas_dnrm2(CD);

    finc = norm_pq_2 * gamma1 * gamma2/(gamma1 + gamma2);

    for (i = 0; i < 10; i++)
        Bx[i] = By[i] = Bz[i] = 0;

    Bxyz(g1->l, g2->l, PA->data[0], PB->data[0], gamma1, \
         g3->l, g4->l, QC->data[0], QD->data[0], gamma2, PQ->data[0], Bx);
    Bxyz(g1->m, g2->m, PA->data[1], PB->data[1], gamma1, \
         g3->m, g4->m, QC->data[1], QD->data[1], gamma2, PQ->data[1], By);
    Bxyz(g1->n, g2->n, PA->data[2], PB->data[2], gamma1, \
         g3->n, g4->n, QC->data[2], QD->data[2], gamma2, PQ->data[2], Bz);

    gsl_vector_free(PA);
    gsl_vector_free(PB);
    gsl_vector_free(QC);
    gsl_vector_free(QD);
    gsl_vector_free(AB);
    gsl_vector_free(CD);
    gsl_vector_free(PQ);

    for (i = 0; i <= L; i++)
        for (j = 0; j <= M; j++)
            for (k = 0; k <= N; k++)
                result += Bx[i] * By[j] * Bz[k] * F_inc_gamma(i+j+k, finc);

    result *= omega(alpha1, alpha2, alpha3, alpha4, norm_ab, norm_cd) \
                * g1->norm * g2->norm * g3->norm * g4->norm \
                * g1->coeff * g2->coeff * g3->coeff * g4->coeff;

    return result;
}

double int2e_basis(const BASIS* b1, const BASIS* b2,
                                const BASIS* b3, const BASIS* b4)
{
    int debug = 0;
    int i, j, k, l, gaussCount_1, gaussCount_2, gaussCount_3, gaussCount_4;
    double result = 0;

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;
    gaussCount_3 = b3->gaussCount;
    gaussCount_4 = b4->gaussCount;

    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            for (k = 0; k < gaussCount_3; k++) {
                for (l = 0; l < gaussCount_4; l++) {
                    result += int2e_gto(&b1->gaussian[i], b1->xyz,
                                                     &b2->gaussian[j], b2->xyz,
                                                     &b3->gaussian[k], b3->xyz,
                                                     &b4->gaussian[l], b4->xyz,
                                                     debug);
                }
            }
        }
    }

    return result;
}

#ifdef __INTEGRAL__INT2E__ONE__
double* int2e_matrix(INPUT_INFO* b)
#else
double**** int2e_matrix(INPUT_INFO* b)
#endif
{
    int i, j, k, l, basis_count, atomCount;
    BASIS *basisSet;

    //INPUT_INFO *b = parse_input(file_name);    

    //ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    atomCount = b->atomCount;

#ifdef __INTEGRAL__INT2E__ONE__
    // use one dimension array output int2e
    double *matrix;
    int I = 0;
    matrix = malloc(sizeof(double) * pow(basis_count, 4));
    for (i = 0; i < basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            for (k = 0; k < basis_count; k++) {
                for (l = 0; l < basis_count; l++) {
                    I = i * pow(basis_count, 3) + j * pow(basis_count, 2) + \
                                                            k * basis_count + l;
                    matrix[I] = int2e_basis(&basisSet[i], 
                                            &basisSet[j], 
                                            &basisSet[k], 
                                            &basisSet[l]);
                }
            }
        }
    }
    return matrix;
#else
    // use four dimension array output int2e
    double ****e2;
    e2 = (double****)malloc(sizeof(double***)*basis_count);
    for (i = 0; i < basis_count; i++) {
        *(e2+i) = (double***)malloc(sizeof(double**)*basis_count);
        for (j = 0; j < basis_count; j++) {
            *(e2[i]+j) = (double**)malloc(sizeof(double*)*basis_count);
            for (k = 0; k < basis_count; k++) {
                *(e2[i][j]+k) = (double*)malloc(sizeof(double)*basis_count);
                for (l = 0; l < basis_count; l++) {
                    e2[i][j][k][l] = int2e_basis(&basisSet[i], 
                                                 &basisSet[j], 
                                                 &basisSet[k], 
                                                 &basisSet[l]);
                }
            }
        }
    }
    return e2;
#endif
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
                    I = i * pow(n, 3) + j * pow(n, 2) + \
                                                            k * n + l;
                    printf("%15.9lf", e[I]);
#else
                    printf("%15.9lf", e[i][j][k][l]);
#endif
                }
            }
            printf("\n");
        }
    }
}