#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "overlap.h"
#include "ints.h"

double fact_l_lambda(int l, int lambda)
{
    // 量子化学中册 P62 (10.6.5)
    return factorial(l) / (factorial(lambda) * factorial(l - lambda));
}

double fi_l_ll_pax_pbx(int ii, int l1, int l2, double pax, double pbx, int flags)
{
    int i, j;
    double sum = 0;

// formula come from 
// Justin T. Fermann and Edward F. Valeev; 
// Fundamentals of Molecular Integrals Evaluation (2.46)

    int top, down;
    top = GSL_MIN(ii, 2*l1-ii);
    down = GSL_MAX(-ii, ii - 2*l2);

    for (i = down; i <= top; i+=2) {
        sum += (fact_l_lambda(l1, (i + ii)/2) * fact_l_lambda(l2, (ii-i)/2) * \
                gsl_pow_int(pax, l1 - (i + ii)/2) * gsl_pow_int(pbx, l2 - (ii-i)/2));
    }

/*
// formula come from PyQuant cints.c
  for (i=0; i<ii+1; i++)
    if ((ii-l1 <= i) && (i <= l2)) 
      sum += fact_l_lambda(l1, ii-i)*fact_l_lambda(l2, i)*pow(pax, l1-ii+i)*pow(pbx,l2-i);

*/

/*
// 《量子化学》中册 P63 第一个公式
    for (i = 0; i <= l1; i++) {
        for (j = 0; j <= l2; j++) {
            if (i + j == ii) {
                sum += (fact_l_lambda(l1, i) * fact_l_lambda(l2, j) *
                                    pow(pax, l1 - i) * pow(pbx, l2 - j));
            }
        }
    }
*/
    return  sum;
}

double I_xyz(int l1, double pax, int l2, double pbx, double gamma, int flags)
{
// 《量子化学》中册 P63 (10.6.8) (10.6.9) (10.6.10)
    int i;
    double sum = 0;

// BUG 2i 
    for (i = 0; i < floor((l1 + l2) * 0.5) + 1; i++) {
        sum += factorial_2(2*i - 1) / gsl_pow_int(2 * gamma, i) * \
                fi_l_ll_pax_pbx(2*i, l1, l2, pax, pbx, flags);
    }
    return sum;
}

double gauss_K(double a, const gsl_vector *A, double b, const gsl_vector *B)
{
// 归一化系数 《量子化学》中册 P58 (10.4.1b)
// A, B 为坐标

    double result, norm_2;
    gsl_vector* v = gsl_vector_alloc(3);

    gsl_vector_memcpy(v, A);
    gsl_vector_sub(v, B);
    norm_2 = gsl_pow_2(gsl_blas_dnrm2(v));

    result = exp(-a * b * norm_2 / (a + b));
    return result;
}

gsl_vector* gaussian_product_center(const double a, const gsl_vector *A, 
                            const double b, const gsl_vector *B, int flags)
{
// Gaussian函数乘积定理计算双中心
// 关于此部分不甚明白
    int i;
    double gamma = a + b;
    //double x1, x2, tmp;

    gsl_vector *center = gsl_vector_alloc(3);
    
    for (i = 0; i < 3; i++)
        center->data[i] = (a * A->data[i] + b * B->data[i]) / gamma;

    // FOR DEBUG
    if (flags == 1) {
        gsl_vector * test = gsl_vector_alloc(3);
        gsl_vector * test2 = gsl_vector_alloc(3);
        gsl_vector_set(test, 0, 0);
        gsl_vector_set(test, 1, 0);
        gsl_vector_set(test, 2, 2.175);
        gsl_vector_memcpy(test2, test);
        gsl_vector_sub(test, B);
        gsl_vector_sub(test2, A);
        if (gsl_vector_max(test) < 1.0E-10 || gsl_vector_max(test2) < 1.0E-10) {
        //printf("----------------------------------------\n");
        vector_output(A, 3, "用于计算中心的第一个坐标:");
        vector_output(B, 3, "用于计算中心的第二个坐标:");
        printf("alpha1 =%10.6lf\talpha2 =%10.6lf\n", a, b);
        vector_output(center, 3, "中心为:");
        }
        gsl_vector_free(test);
        gsl_vector_free(test2);
    }

    return center;
}

double overlap_gto_c(const GTO* g1, const gsl_vector* A, const GTO* g2, const gsl_vector* B, int debug)
{
    double K, gamma, Ix, Iy, Iz, result = 0;
    double normal1, normal2;
    double coeff1, coeff2;
    gsl_vector *PA, *PB, *P;

    PA = gsl_vector_alloc(3);
    PB = gsl_vector_alloc(3);

    normal1 = g1->norm;
    normal2 = g2->norm;
    coeff1 = g1->coeff;
    coeff2 = g2->coeff;

    // compute the coordination of P. 《量子化学》中册, P77
    P = gaussian_product_center(g1->alpha, A, g2->alpha, B, debug);
    gsl_vector_memcpy(PA, P);
    gsl_vector_memcpy(PB, P);

    gsl_vector_sub(PA, A);
    gsl_vector_sub(PB, B);

    gamma = g1->alpha + g2->alpha;

    Ix = I_xyz(g1->l, PA->data[0], g2->l, PB->data[0], gamma, debug);
    Iy = I_xyz(g1->m, PA->data[1], g2->m, PB->data[1], gamma, debug);
    Iz = I_xyz(g1->n, PA->data[2], g2->n, PB->data[2], gamma, debug);

    K = gauss_K(g1->alpha, A, g2->alpha, B);

    result = pow(M_PI/gamma, 1.5) * K * Ix * Iy * Iz * normal1 * normal2 * coeff1 * coeff2;
    if (debug == 2) {
                printf("--------------------------------------------\n");
                vector_output(PA, 3, "PA:");
                vector_output(PB, 3, "PB:");
                gto_output(g1, 1, "g1:");
                gto_output(g2, 1, "g2:");
                printf("%14.8lf%14.8lf%14.8lf%14.8lf%14.8lf\n", 
                                                        K, Ix, Iy, Iz, result);
    }
    gsl_vector_free(P);
    gsl_vector_free(PA);
    gsl_vector_free(PB);
    return result;
}

// 计算两个基函数的重叠积分
double overlap_basis(const BASIS *b1, const gsl_vector *A,
                      const BASIS *b2, const gsl_vector *B, int debug)
{
    //int flags = 0; // debug
    int i, j;
    double result = 0;
    GTO g1, g2;

    // 对组成基组的Gaussian函数循环
    for (i = 0; i < b1->gaussCount; i++) {
        for (j = 0; j < b2->gaussCount; j++) {
            g1 = b1->gaussian[i];
            g2 = b2->gaussian[j];
            result += overlap_gto(&g1, A, &g2, B, debug);
        }
    }
    return result;
}

gsl_matrix* overlap_matrix(INPUT_INFO* b)
{
    int i, j, basis_count;
    BASIS *basisSet;
    double result = 0;
    //double result_check = 0;

    //INPUT_INFO *b = parse_input(file_name);    

    //ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    gsl_matrix *m_overlap = gsl_matrix_calloc(basis_count,basis_count);
    //gsl_matrix *overlap_check = gsl_matrix_calloc(basis_count,basis_count);

    //atom_output((const ATOM_INFO **)alist, b->atomCount);
    //basis_set_output(b->basisSet, b->basisCount, "BASIS");

    for (i = 0; i <  basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            result = overlap_basis(&basisSet[i], basisSet[i].xyz, &basisSet[j],
                            basisSet[j].xyz, 0);
            //result_check = check_overlap(&basisSet[i], &basisSet[j], 0);
            // 设定一个阀值，如果积分值小于某个数就舍去
            if (fabs(result) < 1.0E-12) {
                result = 0;
            }
            gsl_matrix_set(m_overlap, i, j, result);
            //gsl_matrix_set(overlap_check, i, j, result);
        }
    }
    //matrix_output(overlap_check, basis_count, "CHECK OVERLAP:");

    //matrix_output(m_overlap, basis_count, "OVERLAP INTEGRALS:");
    return m_overlap;
}
/*
double check_overlap(const BASIS* b1, const BASIS* b2, int debug)
{
    int i, j;
    int l1, m1, n1, l2, m2, n2;
    double alpha1, alpha2;
    double xa, ya, za, xb, yb, zb;
    int gaussCount_1, gaussCount_2;
    GTO *g1, *g2;
    gsl_vector *A, *B;
    double result = 0;

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;

    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            g1= &b1->gaussian[i];
            g2= &b2->gaussian[j];

            A = b1->xyz;
            B = b2->xyz;

            alpha1 = g1->alpha;
            alpha2 = g2->alpha;


            l1 = g1->l; m1 = g1->m; n1 = g1->n;
            l2 = g2->l; m2 = g2->m; n2 = g2->n;
            xa = gsl_vector_get(A, 0);
            ya = gsl_vector_get(A, 1);
            za = gsl_vector_get(A, 2);

            xb = gsl_vector_get(B, 0);
            yb = gsl_vector_get(B, 1);
            zb = gsl_vector_get(B, 2);

            result += overlap(alpha1, l1, m1, n1, xa, ya, za, 
                              alpha2, l2, m2, n2, xb, yb, zb) * \
                              g1->coeff * g1->norm * g2->coeff * g2->norm;
        }
    }
    return result;
}
*/
