#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_blas.h>
//#include "basis.h"
#include "common.h"
#include "overlap.h"

//#define DEBUG_OVERLAP_SINGLE
//#define DEBUG_DOUBLE_CENTER
//#define DEBUG_OVERLAP_SINGLE

int factorial(int n)
{
    int i, result = 1;

    if (n <= 1) return 1;

    for (i = 1; i <= n; i++)
        result *= i;

    return result;
}

int factorial_2(int n)
{
    int i, result = 1;

    if (n <= 1) return 1;

    for (i = 1; i <= n; i += 2)
        result *= i;

    return result;
}

double fact_l_lambda(int l, int lambda)
{
    // 量子化学中册 P62 (10.6.5)
    return (double)factorial(l) / factorial(lambda) / factorial(l - lambda);
}

double fi_l_ll_pax_pbx(int ii, int l, int ll, double pax, double pbx)
{
// 公式可能有问题，见《量子化学》中册P63第一个公式
    int i, j;
    double sum = 0;

    for (i = 0; i <= l; i++) {
        for (j = 0; j <= ll; j++) {
            if ((i + j) == ii)
                sum += (fact_l_lambda(l, i) * fact_l_lambda(ll, j) * pow(pax, l - i) * pow(pbx, ll - j));        
        }
    }
    return  sum;
}

double I_xyz(double a, int l, double pax, double b, int ll, double pbx)
{
// 《量子化学》中册 P63 (10.6.8) (10.6.9) (10.6.10)
    int i;
    double sum = 0;

    for (i = 0; i <= (l + ll) / 2; i++)
        sum += fi_l_ll_pax_pbx(2*i, l, ll, pax, pbx) * factorial_2(2i - 1) / pow(2 * (a + b), i);
    return sum * sqrt(M_PI / (a + b));
}

double gauss_K(double a, const gsl_vector *A, double b, const gsl_vector *B)
{
// 归一化系数 《量子化学》中册 P58 (10.4.1b)
// A, B 为坐标
    /*
    gsl_vector *AB = gsl_vector_alloc(3);

    gsl_vector_memcpy(AB, A);
    gsl_vector_sub(AB, B);

    double result = exp(-a * b / (a + b) * pow(gsl_blas_dnrm2(AB), 2));
    gsl_vector_free(AB);
    */
    double AB_2 = 0, tmp = 0;
    int i;
    for (i = 0; i < 3; i++) {
        tmp = A->data[i] - B->data[i];
        AB_2 += (tmp * tmp);
    }
    double result = exp((-a) * b / (a + b) * AB_2);
    return result;
}

gsl_vector* gaussian_product_center(const double a, const gsl_vector *A, const double b, const gsl_vector *B)
{
// Gaussian函数乘积定理计算双中心
// 关于此部分不甚明白
    int i;
    double gamma = a + b;
    gsl_vector *center = gsl_vector_alloc(3);
    
    for (i = 0; i < 3; i++)
        center->data[i] = (a * A->data[i] + b * B->data[i]) / gamma;
//#undef DEBUG_GAUSSIAN_PRODUCT
#ifdef DEBUG_GAUSSIAN_PRODUCT
    printf("----------------------------------------\n");
    vector_output(A, 3, "用于计算中心的第一个坐标:");
    vector_output(B, 3, "用于计算中心的第二个坐标:");
    printf("alpha1 =%10.6lf\talpha1 =%10.6lf\n", a, b);
    vector_output(center, 3, "中心为:");
#endif
    return center;
}

double normalize_coeff(const GTO *g)
{
    double alpha = g->alpha;
    double l = g->l;
    double m = g->m;
    double n = g->n;

    return pow(2 * alpha / M_PI, 0.75) * sqrt(pow(4*alpha, l + m + n) / \
        (factorial_2(2*l-1) * factorial_2(2*m-1) * factorial_2(2*n-1)));
}

// 计算两个Gaussian函数的重叠积分
double overlap_single(const GTO *g1, const gsl_vector *A, const GTO *g2, const gsl_vector *B)
{
    double K, gamma, Ix, Iy, Iz, result = 0;
    double normal1, normal2;
    gsl_vector * A_weight = gsl_vector_alloc(3);
    gsl_vector * B_weight = gsl_vector_alloc(3);
    gsl_vector * P = gaussian_product_center(g1->alpha, A, g2->alpha, B);
#ifdef DEBUG_DOUBLE_CENTER
    vector_output(gauss_double_center, 3, "Gaussian双电子中心:\n");
#endif
    gamma = g1->alpha + g2->alpha;
    K = gauss_K(g1->alpha, A, g2->alpha, B);

    // 计算pax, pay, paz, pbx,...
    // 关于此部分不甚明白, 参照PyQuante pyints.py中计算重叠积分部分
    gsl_vector_memcpy(A_weight, A);
    gsl_vector_sub(A_weight, P);
    gsl_vector_memcpy(B_weight, B);
    gsl_vector_sub(B_weight, P);
#ifdef DEBUG_GAUSSIAN_PRODUCT
    vector_output(A_weight, 3, "PA_x_y_z");
    vector_output(B_weight, 3, "PB_x_y_z");
#endif
    Ix = I_xyz(g1->alpha, g1->l, -A_weight->data[0], g2->alpha, g2->l, -B_weight->data[0]);
    Iy = I_xyz(g1->alpha, g1->m, -A_weight->data[1], g2->alpha, g2->m, -B_weight->data[1]);
    Iz = I_xyz(g1->alpha, g1->n, -A_weight->data[2], g2->alpha, g2->n, -B_weight->data[2]);

    normal1 = normalize_coeff(g1);
    normal2 = normalize_coeff(g2);
//#define DEBUG_OVERLAP_SINGLE
#ifdef DEBUG_OVERLAP_SINGLE
    /*printf("K= %8.6lf Ix=%9.6lf Iy=%9.6lf Iy=%9.6lf ", K, Ix, Iy, Iz);
    //printf("N1=%10.6lf N2=%10.6lf Coef1=%10.6lf Coef2=%10.6lf", \
    //    normalize_coeff(g1), normalize_coeff(g2), g1->coeff, g2->coeff);
    */
    //printf("%8.6lf %9.6lf %9.6lf %9.6lf ", K, Ix, Iy, Iz);
    printf("%10.6lf %10.6lf %10.6lf %10.6lf", \
         g1->coeff, g2->coeff, normal1, normal2);
#endif
    result = K * Ix * Iy * Iz * g1->coeff * g2->coeff * normal1 * normal2;
#ifdef DEBUG_OVERLAP_SINGLE
    printf("%10.6lf\n", result);
#endif
    gsl_vector_free(A_weight);
    gsl_vector_free(B_weight);
    return result;
}
/*
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
*/
