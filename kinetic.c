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

double kinetic_single(const GTO *g1, const gsl_vector* A, const GTO* g2, gsl_vector* B, int debug)
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

int main(int argc, char** argv)
{
    int i, j, ii,jj, gauss_i, gauss_j;
    int debug = 0;
    int basis_count_i, basis_count_j;
    ATOM_INFO *atom_list;
    BASIS *basis_i, *basis_j;
    int atom_count;
    double result = 0, result_check = 0;
    gsl_matrix *m_overlap = gsl_matrix_calloc(8,8);
    gsl_matrix *m_overlap_c = gsl_matrix_calloc(8,8);
    int ri = 0, rj = 0;

    char* file_name = "input_file";

    INPUT_INFO *b = parse_input(file_name);    
    //atom_output(b->c, 4);

    atom_list = b->c;
    atom_count = b->n;

    // atom
    for (i = 0; i < atom_count; i++) {
        basis_i = atom_list[i].basis;
        basis_count_i = atom_list[i].basis_count;
        // basis set of atom <i>
        for (ii = 0; ii < basis_count_i; ii++) {
            // atom
            for (j = 0; j < atom_count; j++) {
                basis_j = atom_list[j].basis;
                basis_count_j = atom_list[j].basis_count;
                // basis set of atom <j>
                for (jj = 0; jj < basis_count_j; jj++) {
                    debug = 0;
                    if ((rj == 7 && (ri == 2 || ri == 3 || ri == 4)) || (ri == 7 && (rj == 2 || rj == 3 || rj == 4))) {
                        printf("--------i = %d-----j = %d---------\n", ri,rj);
                        debug = 1;
                    }
                    rj = 0;
                    // gaussian function in basis
                    for (gauss_i = 0; gauss_i < 3; gauss_i++) {
                    // the gaussian function of basis
                        for (gauss_j = 0; gauss_j < 3; gauss_j++) {
                    result += kinetic_single(&basis_i[ii].gaussian[gauss_i], atom_list[i].c, \
                                             &basis_j[jj].gaussian[gauss_j], atom_list[j].c, debug);
               result_check += check_kinetic(&basis_i[ii].gaussian[gauss_i], atom_list[i].c, 
                                             &basis_j[jj].gaussian[gauss_j], atom_list[j].c, 0);
                        }
                    // 设定一个阀值，如果积分值小于某个数就舍去
                    }
                    if (fabs(result) < 1.0E-10)
                        result = 0;
                    gsl_matrix_set(m_overlap, ii, jj, result);
                    gsl_matrix_set(m_overlap_c, ii, jj, result);
                } // end basis of atom <j> integral
            }
        }
    }
    matrix_output(m_overlap, 8, "OVERLAP INTEGRALS:");
    matrix_output(m_overlap_c, 8, "CHECK OVERLAP INTEGRALS:");
    return 0;
}
