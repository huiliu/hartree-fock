#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "basis.h"
#include "common.h"
#include "overlap.h"

int overlap_1(char * basis_base)
{
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


//#define DEBUG_CHECK_INPUT_BASIS
#ifdef DEBUG_CHECK_INPUT_BASIS
    basis_set_output(&b[0], 3, "N-1S:");
    basis_set_output(&b[1], 3, "N-2S:");
    basis_set_output(&b[2], 3, "N-2Px:");
    basis_set_output(&b[3], 3, "N-2Py:");
    basis_set_output(&b[4], 3, "N-2Pz:");
    basis_set_output(&b[5], 3, "H-1S:");
    basis_set_output(&b[6], 3, "H-1S:");
    basis_set_output(&b[7], 3, "H-1S:");
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
                    tmp = overlap_single(&b[u].gaussian[i], &H[0], &b[v].gaussian[j], &H[1]);
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

    matrix_output(m_overlap, basis_count, "重叠积分\n");
    gsl_matrix_free(m_overlap);
#ifdef DEBUG_GAUSS_BASIS_INT
    gsl_matrix_free(gauss_int);
#endif
    gsl_vector_free(h1);
    gsl_vector_free(h2);

    return 0;
}

void overlap_2(char* file_name)
{
    int i, j, ii,jj;
    int gauss_i, gauss_j;
    int basis_count_i, basis_count_j;
    ATOM_INFO *atom_list;
    BASIS *basis_i, *basis_j;
    int atom_count;
    double result, tmp;
    gsl_matrix *m_overlap = gsl_matrix_calloc(8,8);
    int ri = 0, rj = 0;

    INPUT_INFO *b = parse_input(file_name);    
    //atom_output(b->c, 4);

    atom_list = b->c;
    atom_count = b->n;

    for (i = 0; i < atom_count; i++) {
        basis_i = atom_list[i].basis;
        basis_count_i = atom_list[i].basis_count;
        for (ii = 0; ii < basis_count_i; ii++) {
            for (j = 0; j < atom_count; j++) {
                basis_j = atom_list[j].basis;
                basis_count_j = atom_list[j].basis_count;
                for (jj = 0; jj < basis_count_j; jj++) {
                    result = 0;
                    for (gauss_i = 0; gauss_i < 3; gauss_i++) {
                        for (gauss_j = 0; gauss_j < 3; gauss_j++) {
                            tmp = overlap_single(&basis_i[ii].gaussian[gauss_i], atom_list[i].c, &basis_j[jj].gaussian[gauss_j], atom_list[j].c);
                            result += tmp;
                        }
                    }
                    //printf("%15.8lf", result);
                    // 设定一个阀值，如果积分值小于某个数就舍去
                    if (fabs(result) < 1.0E-10)
                        result = 0;
                    gsl_matrix_set(m_overlap, ri, rj, result);
                    if (rj < 7)
                        rj++;
                    else{
                        rj = 0;
                        ri++;
                    }
                }
            }
        }
    }
    matrix_output(m_overlap, 8, "Overlap");

}

int main(int argc, char** argv)
{
    char *basis_base = NULL;
    if (argc < 2)
        basis_base = "basis_set";
    else
        basis_base = argv[1];
    printf("---------Read data from \"%s\"--------------------\n", basis_base);
    overlap_1(basis_base);

    basis_base = "input_file";
    printf("---------Read data from \"%s\"--------------------\n", basis_base);
    overlap_2(basis_base);
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
