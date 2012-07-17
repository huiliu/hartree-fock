#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basis.h"
#include "overlap.h"
#include "common.h"
#include "ints.h"

double check_nuclear(const BASIS* b1, const BASIS* b2, ATOM_INFO **atomList, int atomCount);

double* A_iru(int l1, int l2, double Ax, double Bx, double Cx, double gamma, int debug)
{
// formula (2.18)
    int ii = l1 + l2;
    int i, r, u, I;
    int flags = 0;
    double tmp, *A;

    A = calloc(sizeof(double), ii);

    for (i = 0; i <= ii; i++) {
        for (r = 0; r <= i/2; r++) {
            for (u = 0; u <= (i-2*r)/2; u++) {
                I = i - 2*r -u;
                tmp = pow(-1, u) * factorial(i) * pow(Cx, i - 2*r - 2*u) * \
                        pow(0.25/gamma, r+u) / factorial(r) / factorial(u) / \
                        factorial(i - 2*r - 2*u);
                A[I] += pow(-1, i) * fi_l_ll_pax_pbx(i, l1, l2, Ax, Bx, flags) * tmp;
                if (debug == 2)
                    printf("A_iru = %lf\n", A[I]);
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
                break;
            result += tmp;
        }
        return result * factorial_2(2 * m -1) * exp(-w);;
    }else
        result = factorial_2(2*m -1) / pow(2*w, m + 0.5) * sqrt(M_PI_2);
    return result;
}

double nuclear_elect_attraction_gto(const GTO* g1, const gsl_vector* A, \
        const GTO* g2, const gsl_vector* B, const gsl_vector* C, int debug)
{
    int flags = 0;
    gsl_vector* P = gaussian_product_center(g1->alpha, A, g2->alpha, B, flags);
    gsl_vector *PA, *PB, *PC;
    double gamma;
    double *Ax, *Ay, *Az;
    double norm_pc_2, K;
    double sum = 0, result = 0;
    int i, j, k;
    int l1, n1, m1, l2, m2, n2;

    l1 = g1->l; m1 = g1->m; n1 = g1->n;
    l2 = g2->l; m2 = g2->m; n2 = g2->n;

    gamma = g1->alpha + g2->alpha;
    K = gauss_K(g1->alpha, A, g2->alpha, B);

    PA = gsl_vector_alloc(3);
    PB = gsl_vector_alloc(3);
    PC = gsl_vector_alloc(3);

    gsl_vector_memcpy(PA, P);
    gsl_vector_memcpy(PB, P);
    gsl_vector_memcpy(PC, P);

    gsl_vector_sub(PA, A);
    gsl_vector_sub(PB, B);
    gsl_vector_sub(PC, C);

    norm_pc_2 = pow(gsl_blas_dnrm2(PC), 2);

    Ax = A_iru(l1, l2, PA->data[0], PB->data[0], PC->data[0], gamma, debug);
    Ay = A_iru(m1, m2, PA->data[1], PB->data[1], PC->data[1], gamma, debug);
    Az = A_iru(n1, n2, PA->data[2], PB->data[2], PC->data[2], gamma, debug);

    for (i = 0; i <= l1 + l2; i++)
        for (j = 0; j <= m1 + m2; j++)
            for (k = 0; k <= n1 + n2; k++)
                // 公式与 (2.17)不一样
                sum += Ax[i] * Ay[j] * Az[k] * F_inc_gamma(i+j+k, 
                                                    norm_pc_2*gamma);
    result = 2*M_PI / gamma * K * sum * \
                                g1->norm * g2->norm * g1->coeff * g2->coeff;

if (debug == 2) {
    vector_output(PA, 3, "PA");
    vector_output(PB, 3, "PB");
    vector_output(PC, 3, "PC");
    printf("norm_PC_2 = %20.10lf sum = %20.10lf result = %20.10lf\n", norm_pc_2, sum, result);
}

    gsl_vector_free(PA);
    gsl_vector_free(PB);
    gsl_vector_free(PC);

    return result;
}

// compute the nuclear-electrics attraction integrals of one basis
double nuclear_elect_attraction_basis(const BASIS* b1, const BASIS* b2, 
                                ATOM_INFO **atomList, int atomCount, int debug)
{
    int i, j ,s;
    int gaussCount_1, gaussCount_2;
    ATOM_INFO *atom;
    double result = 0;

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;

if (debug == 2) {
    printf("-----------------------\n");
    basis_set_output(b1, 1, "b1:");
    basis_set_output(b2, 1, "b2:");
}
    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            for (s = 0; s < atomCount; s++) {
                atom = atomList[s];
                result += nuclear_elect_attraction_gto(
                    &b1->gaussian[i], b1->xyz, 
                    &b2->gaussian[j], b2->xyz, 
                    atom->coordination, debug) * atom->n;
            }
if (debug == 2) {
    printf("===============gauss i = %d j = %d================\n", i, j);
    gto_output(&b1->gaussian[i], 1, "i basis:");
    gto_output(&b2->gaussian[j], 1, "j basis:");
}
        }
    }
if (debug == 2)
    printf("\n");

    return result;
}
void nuclear_attraction_matrix(char* file_name)
{
    int debug = 0;
    int i, j, basis_count, atomCount;
    BASIS *basisSet;
    double result = 0, result_check = 0;

    INPUT_INFO *b = parse_input(file_name);    

    ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    atomCount = b->atomCount;

    gsl_matrix *m_overlap = gsl_matrix_calloc(basis_count, basis_count);
    gsl_matrix *m_overlap_c = gsl_matrix_calloc(basis_count, basis_count);

    //atom_output((const ATOM_INFO **)alist, b->atomCount);
    //basis_set_output(b->basisSet, b->basisCount, "BASIS");

    for (i = 0; i <  basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            debug = 0;
            if ( (i == j) && j == 3)
                debug = 2;
            result = nuclear_elect_attraction_basis(&basisSet[i], &basisSet[j],
                            alist, atomCount, debug);
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
    matrix_output(m_overlap, basis_count, "NUCLEAR");
    matrix_output(m_overlap_c, basis_count, "CHECH NUCLEAR");
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
                C = atomList[s]->coordination;

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
