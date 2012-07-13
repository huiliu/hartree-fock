#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "overlap.h"
#include "ints.h"

double kinetic_I_xyz(const GTO* g1, const gsl_vector* A, const GTO* g2, gsl_vector* B, int flags, int debug)
{
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
    kinetic_xyz -= 2 * pow(b, 2) * overlap_gauss(*g1, A, g2_2, B, debug);
    kinetic_xyz -= 0.5 * l2 * (l2 - 1) * overlap_gauss(*g1, A, g2_1, B, debug);

    return kinetic_xyz;
}

// compute the kinetic integral between gaussian function g1 and g2
double kinetic_gto(const GTO *g1, const gsl_vector* A, const GTO* g2, gsl_vector* B, int debug)
{
    double kinetic_I_x= 0;
    double kinetic_I_y= 0;
    double kinetic_I_z= 0;
    double result = 0;

    kinetic_I_x = kinetic_I_xyz(g1, A, g2, B, 0, debug);
    kinetic_I_y = kinetic_I_xyz(g1, A, g2, B, 1, debug);
    kinetic_I_z = kinetic_I_xyz(g1, A, g2, B, 2, debug);

    result = kinetic_I_x + kinetic_I_y + kinetic_I_z;
    return result;
}


double check_kinetic(const GTO *g1, const gsl_vector *A, \
                      const GTO *g2, const gsl_vector *B, int debug)
{
    double alpha1, alpha2, xa, ya, za, xb, yb, zb;
    int l1, m1, n1, l2, m2, n2;
    double result = 0;
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

            result += kinetic(alpha1, l1, m1, n1, xa, ya, za, alpha2, l2, m2, n2, xb, yb, zb);
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
            result += kinetic_gto(&b1->gaussian[i], b1->xyz, &b2->gaussian[i], b2->xyz, debug);
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

    gsl_matrix *m_overlap = gsl_matrix_calloc(8,8);
    gsl_matrix *m_overlap_c = gsl_matrix_calloc(8,8);

    INPUT_INFO *b = parse_input(file_name);    

    //ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    atomCount = b->atomCount;

    for (i = 0; i <  basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            result = kinetic_basis(&basisSet[i], &basisSet[j], debug);
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

int main(int argc, char** argv)
{
    char* file_name = "input_file";

    kinetic_Int_Matrix(file_name);
    return 0;
}
