#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "overlap.h"
#include "hamiltonian.h"
#include "ints.h"


double kinetic_I_xyz(const GTO* g1, const gsl_vector* A, 
                       const GTO* g2, const gsl_vector* B, int flags, int debug)
{
// 徐光宪 《量子化学》 中册 P66 formula (10.7.6)
    double kinetic_xyz = 0;
    double alpha2 = g2->alpha;
    double l2 = 0;

    GTO g2_1 = *g2;
    GTO g2_2 = *g2;

    if (flags == 0) { // 求Px
        l2 = g2->l;
        g2_1.l -= 2;
        g2_2.l += 2;
    }else if (flags == 1) { // 求Py
        l2 = g2->m;
        g2_1.m -= 2;
        g2_2.m += 2;
    }else if (flags == 2) { // 求Pz
        l2 = g2->n;
        g2_1.n -= 2;
        g2_2.n += 2;
    }

    kinetic_xyz = alpha2 * (2*l2 + 1) * overlap_gto(g1, A, g2, B, debug);
    if(!(gtoIsNeg(&g2_2)))  // 貌似没有作用
        kinetic_xyz -= 2 * pow(alpha2, 2) * overlap_gto(g1, A, &g2_2, B, debug);
    if(!(gtoIsNeg(&g2_1)))
        kinetic_xyz -= 0.5 * l2 * (l2 - 1) * overlap_gto(g1, A, &g2_1, B, debug);

    return kinetic_xyz;
/*
// Ohata K, Taketa H, Huzinaga S. Phys Soc Japan, 1966, 2:63
// formula (2.14)

    double result = 0;
    double alpha2 = g2->alpha;
    int l2, m2, n2;

    l2 = g2->l; m2 = g2->m; n2 = g2->n;

    GTO g2_1 = *g2;
    GTO g2_2 = *g2;
    GTO g2_3 = *g2;
    GTO g2_4 = *g2;
    GTO g2_5 = *g2;
    GTO g2_6 = *g2;

    g2_1.l += 2;
    g2_2.m += 2;
    g2_3.n += 2;
    g2_4.l -= 2;
    g2_5.m -= 2;
    g2_6.n -= 2;

    result = alpha2 * (2*(l2+m2+n2) + 3) * overlap_gto(g1, A, g2, B, debug);
    result -= 2 * pow(alpha2, 2) * (overlap_gto(g1, A, &g2_1, B, debug) + 
                                    overlap_gto(g1, A, &g2_2, B, debug) + 
                                    overlap_gto(g1, A, &g2_3, B, debug));
    if(!(gtoIsNeg(&g2_4) || gtoIsNeg(&g2_5) || gtoIsNeg(&g2_6)))
    result -= 0.5 * (l2*(l2-1)*overlap_gto(g1, A, &g2_4, B, debug) + 
                     m2*(m2-1)*overlap_gto(g1, A, &g2_5, B, debug) + 
                     n2*(n2-1)*overlap_gto(g1, A, &g2_6, B, debug));
    return result;
*/
}

// equivalent function kinetic_xyz
double kinetic_I_xyz_c(const GTO* g1, const gsl_vector* A, 
                     const GTO* g2, const gsl_vector* B, int flags, int debug)
{
/* Justin T. Fermann and Edward F. Valeev; 
    Fundamentals of Molecular Integrals Evaluation
    formula (4.13)
 */
    double kinetic_xyz = 0;
    int l1, l2;
    double alpha1, alpha2;

    alpha1 = g1->alpha; 
    alpha2 = g2->alpha; 
    
    GTO g1_1 = *g1;
    GTO g1_2 = *g1;
    GTO g2_1 = *g2;
    GTO g2_2 = *g2;

    if (flags == 0) { // 求Px
    l1 = g1->l; l2 = g2->l;
        g1_1.l--; g1_2.l++;
        g2_1.l--; g2_2.l++;
    }else if (flags == 1) { // 求Px
    l1 = g1->m; l2 = g2->m;
        g1_1.m--; g1_2.m++;
        g2_1.m--; g2_2.m++;
    }else if (flags == 2) { // 求Px
    l1 = g1->n; l2 = g2->n;
        g1_1.n--; g1_2.n++;
        g2_1.n--; g2_2.n++;
    }

        kinetic_xyz = 0.5 * l1 * l2 * overlap_gto(&g1_1, A, &g2_1, B, debug);
        kinetic_xyz += 2 * alpha1 * alpha2 * overlap_gto(&g1_2, A, &g2_2, B, debug);
        kinetic_xyz -= alpha1 * l2 * overlap_gto(&g1_2, A, &g2_1, B, debug);
        kinetic_xyz -= alpha2 * l1 * overlap_gto(&g1_1, A, &g2_2, B, debug);

    if (debug == 1) {
        printf("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n");  
        gto_output(&g1_1, 1, "g1_1:");
        gto_output(&g1_2, 1, "g1_2:");
        gto_output(&g2_1, 1, "g2_1:");
        gto_output(&g2_2, 1, "g2_2:");
        printf("          kinetic_xyz = %lf\n", kinetic_xyz);
    }

    return kinetic_xyz;
}

// compute the kinetic integral between gaussian function g1 and g2
double kinetic_gto(const GTO *g1, const gsl_vector* A, 
                   const GTO* g2, const gsl_vector* B, int debug)
{
//  徐光宪 《量子化学》 中册 P66 formula (10.7.4)
    double kinetic_I_x= 0;
    double kinetic_I_y= 0;
    double kinetic_I_z= 0;
    double result = 0;

    if (debug == 2) {
        printf("^^^^^^^ ^^^^^^^^ ^^^^^^^^ ^^^^^^^^^^^ ^^^^^^^^\n");
        vector_output(A, 3, "第一个原子坐标:");
        gto_output(g1, 1, "基函数:");
        vector_output(B, 3, "第二个原子坐标:");
        gto_output(g2, 1, "基函数:");
        printf("Ix%15.9lf%15.9lf%15.9lf\n", kinetic_I_x, kinetic_I_y, kinetic_I_z);
    }

    kinetic_I_x = kinetic_I_xyz(g1, A, g2, B, 0, debug);
    kinetic_I_y = kinetic_I_xyz(g1, A, g2, B, 1, debug);
    kinetic_I_z = kinetic_I_xyz(g1, A, g2, B, 2, debug);

    result = kinetic_I_x + kinetic_I_y + kinetic_I_z;
    return result;
}

// compute the kinetic integral between basis b1 and b2
double kinetic_basis(const BASIS* b1, const BASIS* b2, int debug)
{
    int i, j;
    int gaussCount_1, gaussCount_2;
    double result = 0;

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;

    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            result += kinetic_gto(&b1->gaussian[i], b1->xyz, 
                                  &b2->gaussian[j], b2->xyz, debug);
        }
    }

    return result;
}

gsl_matrix* kinetic_matrix(INPUT_INFO* b)
{
    int debug = 0;
    int i, j, basis_count, atomCount;
    BASIS *basisSet;
    double result = 0;
    //double result_check = 0;

    //INPUT_INFO *b = parse_input(file_name);    

    //ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    atomCount = b->atomCount;

    gsl_matrix *m_kinetic = gsl_matrix_calloc(basis_count, basis_count);
    //gsl_matrix *kinetic_check = gsl_matrix_calloc(basis_count, basis_count);

    for (i = 0; i <  basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            debug = 0;
            result = kinetic_basis(&basisSet[i], &basisSet[j], debug);
            //result_check = check_kinetic(&basisSet[i], &basisSet[j], 0);
            // 设定一个阀值，如果积分值小于某个数就舍去
            if (fabs(result) < 1.0E-12) {
                result = 0;
            }
            gsl_matrix_set(m_kinetic, i, j, result);
            //gsl_matrix_set(kinetic_check, i, j, result_check);
        }
    }
    //matrix_output(m_kinetic, basis_count, "KINETIC");
    //matrix_output(kinetic_check, basis_count, "CHECK KINETIC:");
    return m_kinetic;
}

// --------------------------------------------------------------------------
// COMPUTE NUCLEAR - ELECTRON ATTRACTION INTERACTION
//

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

gsl_matrix* nuclear_attraction_matrix(INPUT_INFO* b)
{
    int debug = 0;
    int i, j, basis_count, atomCount;
    BASIS *basisSet;
    double result = 0;

    //INPUT_INFO *b = parse_input(file_name);    

    ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    atomCount = b->atomCount;

    gsl_matrix *m_attraction = gsl_matrix_calloc(basis_count, basis_count);
    //gsl_matrix *m_overlap_c = gsl_matrix_calloc(basis_count, basis_count);

    //atom_output((const ATOM_INFO **)alist, b->atomCount);
    //basis_set_output(b->basisSet, b->basisCount, "BASIS");

    for (i = 0; i <  basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            /*
            debug = 0;
            if ( (i == j) && j == 3)
                debug = 2;
            */
            result = nuclear_elect_attraction_basis(&basisSet[i], &basisSet[j],
                            alist, atomCount, debug);
            //result_check = check_nuclear(&basisSet[i], &basisSet[j],
            //                alist, atomCount);
            // 设定一个阀值，如果积分值小于某个数就舍去
            if (fabs(result) < 1.0E-12) {
                result = 0;
            }
            gsl_matrix_set(m_attraction, i, j, result);
            //gsl_matrix_set(m_overlap_c, i, j, result_check);
        }
    }
    //matrix_output(m_overlap, basis_count, "NUCLEAR");
    //matrix_output(m_overlap_c, basis_count, "CHECH NUCLEAR");
    return m_attraction;
}

// hamiltonian matrix
gsl_matrix* hamiltonian(INPUT_INFO* b)
{
    int n = b->basisCount;
    gsl_matrix* a = nuclear_attraction_matrix(b);
    gsl_matrix* k = kinetic_matrix(b);
    gsl_matrix* h = gsl_matrix_alloc(n, n);

    gsl_matrix_memcpy(h, k);
    gsl_matrix_sub(h, a);

    gsl_matrix_free(a);
    gsl_matrix_free(k);

    return h;
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

double check_kinetic(const BASIS* b1, const BASIS* b2, int debug)
{
    int i, j;
    int gaussCount_1, gaussCount_2;
    double alpha1, alpha2, xa, ya, za, xb, yb, zb;
    int l1, m1, n1, l2, m2, n2;
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

            result += kinetic(alpha1, l1, m1, n1, xa, ya, za, 
                              alpha2, l2, m2, n2, xb, yb, zb) * \
                              g1->coeff * g1->norm * g2->coeff * g2->norm;
        }
    }
    return result;
}
