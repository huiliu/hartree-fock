#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "overlap.h"
#include "ints.h"

int gtoIsNeg(const GTO* g)
{
    if (g->l < 0 || g->m < 0 || g->n < 0)
        return 1;
    else
        return 0;
}

double kinetic_I_xyz(const GTO* g1, const gsl_vector* A, 
                     const GTO* g2, gsl_vector* B, int flags, int debug)
{
// 徐光宪 《量子化学》 中册 P66 formula (10.7.6)
    double kinetic_xyz = 0;
    double b = g2->alpha;
    double l2 = 0;

    GTO g2_1 = *g2;
    GTO g2_2 = *g2;

    if (flags == 0) { // 求Px
        l2 = g2->l;
        g2_1.l -= 2;
        g2_2.l += 2;
    }else if (flags == 1) { // 求Px
        l2 = g2->m;
        g2_1.m -= 2;
        g2_2.m += 2;
    }else if (flags == 2) { // 求Px
        l2 = g2->n;
        g2_1.n -= 2;
        g2_2.n += 2;
    }

    kinetic_xyz += b * (2*l2 + 1) * overlap_gauss(*g1, A, *g2, B, debug);
    if(!(gtoIsNeg(&g2_2)))  // 貌似没有作用
        kinetic_xyz -= 2 * pow(b, 2) * overlap_gauss(*g1, A, g2_2, B, debug);
    if(!(gtoIsNeg(&g2_1)))
        kinetic_xyz -= 0.5 * l2 * (l2 - 1) * overlap_gauss(*g1, A, g2_1, B, debug);

    return kinetic_xyz;
}

// equivalent function kinetic_xyz
double kinetic_I_xyz_c(const GTO* g1, const gsl_vector* A, 
                     const GTO* g2, gsl_vector* B, int flags, int debug)
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

        kinetic_xyz += 0.5 * l1 * l2 * overlap_gauss(g1_1, A, g2_1, B, debug);
        kinetic_xyz += 2 * alpha1 * alpha2 * overlap_gauss(g1_2, A, g2_2, B, debug);
        kinetic_xyz -= alpha1 * l2 * overlap_gauss(g1_2, A, g2_1, B, debug);
        kinetic_xyz -= alpha2 * l1 * overlap_gauss(g1_1, A, g2_2, B, debug);

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
                   const GTO* g2, gsl_vector* B, int debug)
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
                              alpha2, l2, m2, n2, xb, yb, zb);
        }
    }
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

void kinetic_Int_Matrix(const char* file_name)
{
    int debug = 0;
    int i, j, basis_count, atomCount;
    BASIS *basisSet;
    double result, result_check = 0;

    INPUT_INFO *b = parse_input(file_name);    

    //ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    atomCount = b->atomCount;

    gsl_matrix *m_overlap = gsl_matrix_calloc(basis_count, basis_count);
    gsl_matrix *m_overlap_c = gsl_matrix_calloc(basis_count, basis_count);

    for (i = 0; i <  basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            debug = 0;
            result = kinetic_basis(&basisSet[i], &basisSet[j], debug);
            result_check = check_kinetic(&basisSet[i], &basisSet[j], 0);
            // 设定一个阀值，如果积分值小于某个数就舍去
            if (fabs(result) < 1.0E-12) {
                result = 0;
                result_check = 0;
            }
            gsl_matrix_set(m_overlap, i, j, result);
            gsl_matrix_set(m_overlap_c, i, j, result_check);
        }
    }
    matrix_output(m_overlap, basis_count, "KINETIC");
    matrix_output(m_overlap_c, basis_count, "CHECH KINETIC");
}

int main(int argc, char** argv)
{
    char *basis_base = NULL;

    if (argc < 2)
        basis_base = "input_file";
    else
        basis_base = argv[1];

    kinetic_Int_Matrix(basis_base);
    return 0;
}
