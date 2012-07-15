#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basis.h"
#include "overlap.h"
#include "common.h"
#include "ints.h"

double check_nuclear(const BASIS* b1, const BASIS* b2, ATOM_INFO **atomList, int atomCount);

double* A_iru(int l1, int l2, double Ax, double Bx, double Cx, double gamma)
{
// formula (2.18)
    int ii = l1 + l2;
    int i, r, u, I;
    int flags = 0;
    double tmp, *A;

    A = calloc(sizeof(double), ii);

    for (i = 0; i <= ii; i++) {
        int rr = i / 2;
        for (r = 0; r <= rr; r++) {
            int uu = (i - 2*r) / 2;
            for (u = 0; u <= uu; u++) {
                I = i - 2*r -u;
                tmp = pow(-1, u) * factorial(i) * pow(Cx, i - 2*r - 2*u) * \
                        pow(0.25/gamma, r+u) / factorial(r) / factorial(u) / \
                        factorial(i - 2*r - 2*u);
                A[I] += pow(-1, i) * fi_l_ll_pax_pbx(i, l1, l2, Ax, Bx, flags) * tmp;
            }
        }
    }
    return A;
}

#define F_INC_GAMMA_CYCLE    100
#define F_INC_GAMMA_delta  1.0E-12
double F_inc_gamma(int m ,double w)
{
    double result = 0;
    double tmp = 0;
    int i;
    
    if (w < 17) {
        result = tmp = 1.0 / factorial_2(2*m + 1);
        for (i = 1; i < F_INC_GAMMA_CYCLE; i++) {
            tmp *= ((2*w) / (2*m + 2*i + 1));
            if ((tmp - F_INC_GAMMA_delta) < 0)
                return result * factorial_2(2 * m -1) * exp(-w);;
            result += tmp;
        }
        return result * factorial_2(2 * m -1) * exp(-w);;
    }else
        result = factorial_2(2*m -1) / pow(2*w, m + 0.5) * sqrt(M_PI_2);
    return result;
}

double nuclear_elect_attraction_gto(const GTO* g1, const gsl_vector* A, \
        const GTO* g2, const gsl_vector* B, const gsl_vector* C)
{
    int flags = 0;
    gsl_vector* P = gaussian_product_center(g1->alpha, A, g2->alpha, B, flags);
    gsl_vector *PA, *PB, *PC;
    double gamma;
    double **Axyz;
    double norm_pc_2, K, normlize_coeff_a, normlize_coeff_b;
    double sum = 0, result = 0;
    int i, j, k;
    int l1, n1, m1, l2, m2, n2;

    Axyz = (double **)malloc(sizeof(double *)*3);
    l1 = g1->l; m1 = g1->m; n1 = g1->n;
    l2 = g2->l; m2 = g2->m; n2 = g2->n;

    gamma = g1->alpha + g2->alpha;
    K = gauss_K(g1->alpha, A, g2->alpha, B);
    normlize_coeff_a = normalize_coeff(g1);
    normlize_coeff_b = normalize_coeff(g2);

    PA = gsl_vector_alloc(3);
    PB = gsl_vector_alloc(3);
    PC = gsl_vector_alloc(3);

    gsl_vector_memcpy(PA, P);
    gsl_vector_memcpy(PB, P);
    gsl_vector_memcpy(PC, P);

    norm_pc_2 = pow(gsl_blas_dnrm2(PC), 2);
    //norm_pc_2 = pow(PC->data[0], 2) + pow(PC->data[1], 2) + pow(PC->data[2], 2);

    gsl_vector_sub(PA, A);
    gsl_vector_sub(PB, A);
    gsl_vector_sub(PC, A);


    for (i = 0; i < 3; i++)
        Axyz[i] = A_iru(g1->l, g2->l, PA->data[i], PB->data[i], PC->data[i], gamma);

/*
\* 等同于上部分
    Ax = A_iru(g1->l, g2->l, P->data[0] - A->data[0], P->data[0] - B->data[0], P->data[0] - C->data[0], gamma);
    Ax = A_iru(g1->l, g2->l, P->data[1] - A->data[1], P->data[1] - B->data[1], P->data[1] - C->data[1], gamma);
    Ax = A_iru(g1->l, g2->l, P->data[2] - A->data[2], P->data[2] - B->data[2], P->data[2] - C->data[2], gamma);
*/
    for (i = 0; i <= l1 + l2; i++)
        for (j = 0; j <= m1 + m2; j++)
            for (k = 0; k <= n1 + n2; k++)
                // 公式与 (2.17)不一样
                sum += Axyz[0][i] * Axyz[1][j] * Axyz[2][k] * F_inc_gamma(i+j+k, norm_pc_2*gamma);
    result = 2*M_PI / gamma * K * sum * normlize_coeff_a * normlize_coeff_b;

    gsl_vector_free(PA);
    gsl_vector_free(PB);
    gsl_vector_free(PC);

    return result;
}

// compute the nuclear-electrics attraction integrals of one basis
double nuclear_elect_attraction_basis(const BASIS* b1, const BASIS* b2, ATOM_INFO **atomList, int atomCount)
{
    int i, j ,s;
    int gaussCount_1, gaussCount_2;
    double result = 0;

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;

    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            for (s = 0; s < atomCount; s++) {
                result += nuclear_elect_attraction_gto(&b1->gaussian[i], b1->xyz, &b2->gaussian[j], b2->xyz, atomList[s]->coordination);
            }
        }
    }

    return result;
}
void nuclear_attraction_matrix(char* file_name)
{
    int i, j, basis_count, atomCount;
    BASIS *basisSet;
    double result, result_check = 0;

    gsl_matrix *m_overlap = gsl_matrix_calloc(8,8);
    gsl_matrix *m_overlap_c = gsl_matrix_calloc(8,8);

    INPUT_INFO *b = parse_input(file_name);    

    ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    atomCount = b->atomCount;

    //atom_output((const ATOM_INFO **)alist, b->atomCount);
    //basis_set_output(b->basisSet, b->basisCount, "BASIS");

    for (i = 0; i <  basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            result = nuclear_elect_attraction_basis(&basisSet[i], &basisSet[j],
                            alist, atomCount);
            result_check = check_nuclear(&basisSet[i], &basisSet[j],
                            alist, atomCount);
            // 设定一个阀值，如果积分值小于某个数就舍去
            if (fabs(result) < 1.0E-12) {
                result = 0;
                result_check = 0;
            }
            gsl_matrix_set(m_overlap, i, j, result);
            gsl_matrix_set(m_overlap_c, i, j, result_check);
        }
    }
    matrix_output(m_overlap, 8, "NUCLEAR");
    matrix_output(m_overlap_c, 8, "CHECH NUCLEAR");
}

double check_nuclear(const BASIS* b1, const BASIS* b2, ATOM_INFO **atomList, int atomCount)
{
    int i, j, s;
    int gaussCount_1, gaussCount_2;
    double alpha1, alpha2, xa, ya, za, xb, yb, zb, xc, yc, zc, norm1, norm2;
    int l1, m1, n1, l2, m2, n2;
    gsl_vector *A, *B, *C;
    GTO *g1, *g2;
    double result = 0;

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;

    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            g1 = &b1->gaussian[i];
            g2 = &b2->gaussian[i];

            A = b1->xyz;
            B = b2->xyz;

            l1 = g1->l; m1 = g1->m; n1 = g1->n; alpha1 = g1->alpha; norm1 = g1->norm;
            l2 = g2->l; m2 = g2->m; n2 = g2->n; alpha2 = g2->alpha; norm2 = g2->norm;

            xa = gsl_vector_get(A, 0);
            ya = gsl_vector_get(A, 1);
            za = gsl_vector_get(A, 2);

            xb = gsl_vector_get(B, 0);
            yb = gsl_vector_get(B, 1);
            zb = gsl_vector_get(B, 2);

            for (s = 0; s < atomCount; s++) {
                C = atomList[i]->coordination;

                xc = gsl_vector_get(C, 0);
                yc = gsl_vector_get(C, 1);
                zc = gsl_vector_get(C, 2);

                result = nuclear_attraction(xa, ya, za, norm1, l1, m1, n1, alpha1, \
                                xb, yb, zb, norm2, l2, m2, n2, alpha2, \
                                xc, yc, zc);
            }
        }
    }
    return result;
}
int main(int argc, char** argv)
{
    char *basis_base = NULL;

    if (argc < 2)
        basis_base = "input_file";
    else
        basis_base = argv[1];

    nuclear_attraction_matrix(basis_base);
    return 0;
}
/*
$COORD
N   7   0.30926 0.30926 0.30926
H   1   2.175   0.0     0.0
H   1   0.0     2.175   0.0
H   1   0.0     0.0     2.175
$END
$BASIS
STO-3G
****
N   0
S   3   1.00
     0.9910616896E+02  0.1543289673E+00
     0.1805231239E+02  0.5353281423E+00
     0.4885660238E+01  0.4446345422E+00
SP   3   1.00
     0.3780455879E+01 -0.9996722919E-01  0.1559162750E+00
     0.8784966449E+00  0.3995128261E+00  0.6076837186E+00
     0.2857143744E+00  0.7001154689E+00  0.3919573931E+00
****
H    0
S   3   1.00
     0.3425250914E+01  0.1543289673E+00
     0.6239137298E+00  0.5353281423E+00
     0.1688554040E+00  0.4446345422E+00
****
H   0
S   3   1.00
     0.3425250914E+01  0.1543289673E+00
     0.6239137298E+00  0.5353281423E+00
     0.1688554040E+00  0.4446345422E+00
****
H   0
S   3   1.00
     0.3425250914E+01  0.1543289673E+00
     0.6239137298E+00  0.5353281423E+00
     0.1688554040E+00  0.4446345422E+00
****
$END
*/
