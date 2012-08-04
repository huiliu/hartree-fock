#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "overlap.h"
#include "basis.h"
#include "int2e.h"
#include "eri.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

double theta(int l, int l1, int l2, double PA, double PB, int r, double gamma)
{
// David B. Cook; Handbook of Computational Quantum Chemistry. P249
    int debug = 0;
    return fi_l_ll_pax_pbx(l, l1, l2, PA, PB, debug) * factorial(l) * \
            gsl_pow_int(gamma, r-l) / factorial(r) / factorial(l - 2*r);
}

double B(int l, int l1, int l2, double PA, double PB, int r, double gamma1,
         int ll,int l3, int l4, double QC, double QD, int rr,double gamma2,
         int i, double px)
{
    double result = 0;
    double delta = 0.25 * (1/gamma1 + 1/gamma2);

    result = gsl_pow_int(-1, l+i) * theta(l, l1, l2, PA, PB, r, gamma1) * \
                            theta(ll, l3, l4, QC, QD, rr, gamma2);
    result *= gsl_pow_int(delta, 2*(r+rr)+i-l-ll) * factorial(l + ll -2*(r + rr)) * \
              gsl_pow_int(px, l + ll - 2*(r + rr + i)) * gsl_pow_int(4, r + rr - l - ll);
    result /= (factorial(i) * factorial(l + ll - 2*(r + rr + i)));

    return result;
}

double omega(double alpha1, double alpha2, double alpha3, double alpha4, 
                                                        double AB, double CD)
{
    double gamma1 = alpha1 + alpha2;
    double gamma2 = alpha3 + alpha4;

    return 2*gsl_pow_2(M_PI)/(gamma1*gamma2) * sqrt(M_PI / (gamma1+gamma2)) * \
                exp(-alpha1*alpha2 * gsl_pow_2(AB) / gamma1 
                                       -alpha3*alpha4 * gsl_pow_2(CD) / gamma2);
}

void Bxyz(int l1, int l2, double PA, double PB, double gamma1, 
          int l3, int l4, double QC, double QD, double gamma2,
          double px, double* result)
{
    int l, r, i, ll, rr, I;


    for (l = 0; l <= l1+l2; l++) {
        for (r = 0; r <= l*0.5; r++) {
            for (ll = 0; ll <= l3 + l4; ll++) {
                for (rr = 0; rr <= ll*0.5; rr++) {
                    for (i = 0; i <= (l+ll - 2*(r+rr))*0.5; i++) {
                        I = l + ll - 2*(r + rr) - i;
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
    PQ = gsl_vector_alloc(3);

    gsl_vector_memcpy(PA, P); gsl_vector_memcpy(PB, P);
    gsl_vector_memcpy(QC, Q); gsl_vector_memcpy(QD, Q);
    gsl_vector_memcpy(AB, A); gsl_vector_memcpy(CD, C);
    gsl_vector_memcpy(PQ, P);

    gsl_vector_sub(PA, A); gsl_vector_sub(PB, B);
    gsl_vector_sub(QC, C); gsl_vector_sub(QD, D);
    gsl_vector_sub(AB, B); gsl_vector_sub(CD, D);
    gsl_vector_sub(PQ, Q);


    norm_ab = gsl_blas_dnrm2(AB);
    norm_cd = gsl_blas_dnrm2(CD);

    norm_pq_2 = gsl_pow_2(gsl_blas_dnrm2(PQ));

    finc = norm_pq_2 * gamma1 * gamma2/(gamma1 + gamma2);

    for (i = 0; i < 10; i++)
        Bx[i] = By[i] = Bz[i] = 0;

    Bxyz(g1->l, g2->l, PA->data[0], PB->data[0], gamma1, \
         g3->l, g4->l, QC->data[0], QD->data[0], gamma2, PQ->data[0], Bx);

    Bxyz(g1->m, g2->m, PA->data[1], PB->data[1], gamma1, \
         g3->m, g4->m, QC->data[1], QD->data[1], gamma2, PQ->data[1], By);

    Bxyz(g1->n, g2->n, PA->data[2], PB->data[2], gamma1, \
         g3->n, g4->n, QC->data[2], QD->data[2], gamma2, PQ->data[2], Bz);


    for (i = 0; i <= L; i++) {
        for (j = 0; j <= M; j++) {
            for (k = 0; k <= N; k++) {
                if (debug == 1)
                    printf("%12.6lf%12.6lf%12.6lf %d %d %d\n", Bx[i], By[j], Bz[k], g4->l, g4->m, g4->n);
                result += Bx[i] * By[j] * Bz[k] * F_inc_gamma(i+j+k, finc);
            }
        }
    }

    result *= omega(alpha1, alpha2, alpha3, alpha4, norm_ab, norm_cd) \
                * g1->norm * g2->norm * g3->norm * g4->norm \
                * g1->coeff * g2->coeff * g3->coeff * g4->coeff;

    gsl_vector_free(PA);
    gsl_vector_free(PB);
    gsl_vector_free(QC);
    gsl_vector_free(QD);
    gsl_vector_free(AB);
    gsl_vector_free(CD);
    gsl_vector_free(PQ);
    gsl_vector_free(P);
    gsl_vector_free(Q);

    return result;
}

double int2e_basis(const BASIS* b1, const BASIS* b2,
                   const BASIS* b3, const BASIS* b4, int debug)
{
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
    int debug = 0;
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
    matrix = malloc(sizeof(double) * gsl_pow_4(basis_count));
    for (i = 0; i < basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            for (k = 0; k < basis_count; k++) {
                for (l = 0; l < basis_count; l++) {
                    I = i * gsl_pow_3(basis_count) + j * gsl_pow_2(basis_count) + \
                                                            k * basis_count + l;
                    matrix[I] = int2e_basis(&basisSet[i],
                                            &basisSet[j], 
                                            &basisSet[k], 
                                            &basisSet[l], debug);
                }
            }
        }
    }
    return matrix;
#else
    // use four dimension array output int2e
    double ****e2;
    e2 = (double****)Malloc(sizeof(double***)*basis_count);
    for (i = 0; i < basis_count; i++) {
        *(e2+i) = (double***)Malloc(sizeof(double**)*basis_count);
        for (j = 0; j < basis_count; j++) {
            *(e2[i]+j) = (double**)Malloc(sizeof(double*)*basis_count);
            for (k = 0; k < basis_count; k++) {
                *(e2[i][j]+k) = (double*)Calloc(sizeof(double), basis_count);
            }
        }
    }

    omp_set_num_threads(2);
    //#pragma omp parallel for private(j, k, l)
    for (i = 0; i < basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            for (k = 0; k < basis_count; k++) {
                for (l = 0; l < basis_count; l++) {
                    //debug = 0;
                    //if (i == 0 && j == 2 && k == 0 && l == 2)
                    //    debug = 1;
                    //if (e2[i][j][k][l] != 0)
                    if (chkSYM(e2, i, j, k, l))
                        continue;
                    e2[i][j][k][l] = \
                    e2[i][j][l][k] = \
                    e2[j][i][k][l] = \
                    e2[j][i][l][k] = \
                    e2[k][l][i][j] = \
                    e2[k][l][j][i] = \
                    e2[l][k][i][j] = \
                    e2[l][k][j][i] = int2e_basis(&basisSet[i], 
                                                 &basisSet[j], 
                                                 &basisSet[k], 
                                                 &basisSet[l], debug);
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

// check the symtery of two-electron integral
int chkSYM(double ****e, int i, int j, int k, int l)
{
    //if (fabs(e[i][j][k][l]) > 1.0E-10)
    if (e[i][j][k][l] != 0)
        return 1;
    else if (l < k && e[i][j][l][k] == 0) {e[i][j][k][l] == 0; return 1;}
    else if (j < i && e[j][i][l][k] == 0) {e[i][j][k][l] == 0; return 1;}
    else if (k < i && l < j && e[k][l][i][j] == 0) {e[i][j][k][l] == 0; return 1;}
    else if (l < i && k < j && l <= k && e[l][k][i][j] == 0) {e[i][j][k][l] == 0; return 1;}
    else
        return 0;
}
