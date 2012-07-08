#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "overlap.h"

double kinetic_I_xyz(const GTO* g1, const gsl_vector* A, const GTO* g2, gsl_vector* B, int flags)
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

    kinetic_xyz = -0.5 * l2 * (l2 - 1) * overlap_single(g1, A, &g2_2, B);
    kinetic_xyz += b * (2*l2 + 1) * overlap_single(g1, A, g2, B);
    kinetic_xyz -= 2 * pow(b, 2) * overlap_single(g1, A, &g2_1, B);

    return kinetic_xyz;
}

double kinetic_single(const GTO *g1, const gsl_vector* A, const GTO* g2, gsl_vector* B)
{
    double kinetic_I_x= 0;
    double kinetic_I_y= 0;
    double kinetic_I_z= 0;

    kinetic_I_x = kinetic_I_xyz(g1, A, g2, B, 0);
    kinetic_I_y = kinetic_I_xyz(g1, A, g2, B, 1);
    kinetic_I_z = kinetic_I_xyz(g1, A, g2, B, 2);

    return (kinetic_I_x + kinetic_I_y + kinetic_I_z);
}

int main(int argc, char** argv)
{
    char *basis_base = NULL;
    if (argc < 2)
        basis_base = "basis_set";
    else
        basis_base = argv[1];

    BASIS *b = read_basis(basis_base);    

    double result = 0, tmp;
    int i, j, u, v;
    int basis_count = 8;
    // 重叠积分不位置有关，基组处于各自的中心，
    // 通过记录每个原子所包含的基函数数目来判断某个基函数处于哪里
    int atom_1_basis_num, atom_2_basis_num;
    gsl_matrix *m_overlap = gsl_matrix_alloc(basis_count, basis_count);
    gsl_vector H[2];
    gsl_vector *N = gsl_vector_calloc(3);
    gsl_vector *h1 = gsl_vector_calloc(3);
    gsl_vector *h2 = gsl_vector_calloc(3);
    gsl_vector *h3 = gsl_vector_calloc(3);

    N->data[0] = 0.30926;
    N->data[1] = 0.30926;
    N->data[2] = 0.30926;

    h1->data[0] = 2.17510;
    h1->data[1] = 0.0;
    h1->data[2] = 0.0;

    h2->data[0] = 0.0;
    h2->data[1] = 2.17510;
    h2->data[2] = 0.0;

    h3->data[0] = 0.0;
    h3->data[1] = 0.0;
    h3->data[2] = 2.17510;

#ifdef DEBUG_GAUSS_BASIS_INT
    gsl_matrix *gauss_int = gsl_matrix_calloc(3, 3);
#endif


#undef DEBUG_CHECK_INPUT_BASIS
#ifdef DEBUG_CHECK_INPUT_BASIS
    basis_set_output(&b[0], 3, "N-1S:\n");
    basis_set_output(&b[1], 3, "N-2S:\n");
    basis_set_output(&b[2], 3, "N-2P:\n");
    basis_set_output(&b[3], 3, "N-2P:\n");
    basis_set_output(&b[4], 3, "N-2P:\n");
    basis_set_output(&b[5], 3, "H-1S:\n");
    basis_set_output(&b[6], 3, "H-1S:\n");
    basis_set_output(&b[7], 3, "H-1S:\n");
#endif

    atom_1_basis_num = 5;
    atom_2_basis_num = 3;   //暂时只考虑N原子的一个P轨道
    // 基函数的数目
    for (u = 0; u < basis_count; u++) {
        for (v = 0; v < basis_count; v++) {
            result = 0;

            // 判断某个基组属于哪个中心
            // 此处基函数是按顺序存放在一维数组中
            if (u < atom_1_basis_num) {
            // 0 - 4为N的基函数
                H[0] = *N;
            }else if (u == 5){
                H[0] = *h1;
            }else if (u == 6){
                H[0] = *h2;
            }else if (u == 7){
                H[0] = *h3;
            }

            if (v < atom_1_basis_num) {
            // 0 - 4为N的基函数
                H[1] = *N;
            }else if (v == 5){
                H[1] = *h1;
            }else if (v == 6){
                H[1] = *h2;
            }else if (v == 7){
                H[1] = *h3;
            }

#undef DEBUG_BASIS_COOR
#ifdef DEBUG_BASIS_COOR
    printf("-------------------%d-%d---------------------------\n", u, v);
            basis_set_output(&b[u], 3, "第一个基函数参数:\n");
            basis_set_output(&b[v], 3, "第二个基函数参数:\n");
            vector_output(&H[0], 3, "第一个H原子的坐标为:\n");
            vector_output(&H[1], 3, "第二个H原子的坐标为:\n");
#endif
            // Gaussian函数的数目为3
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    tmp = kinetic_single(&b[u].gaussian[i], &H[0], &b[v].gaussian[j], &H[1]);
                    result += tmp;
#ifdef DEBUG_GAUSS_BASIS_INT
                    gsl_matrix_set(gauss_int, i, j, tmp);
#endif
                }
            }
#ifdef DEBUG_GAUSS_BASIS_INT
            matrix_output(gauss_int, 3, "Gaussian函数的积分\n");
#endif
            // 设定一个阀值，如果积分值小于某个数就舍去
            if (fabs(result) < 1.0E-10)
                result = 0;

            gsl_matrix_set(m_overlap, u, v, result);
//            printf("u = %d\tv = %d%20.10lf\n", u, v, result);
        }
    }

    matrix_output(m_overlap, basis_count, "动能积分\n");
    gsl_matrix_free(m_overlap);
#ifdef DEBUG_GAUSS_BASIS_INT
    gsl_matrix_free(gauss_int);
#endif
    gsl_vector_free(h1);
    gsl_vector_free(h2);

    return 0;
}
