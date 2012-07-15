#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

    if (n % 2 == 0) {
        if (n <= 2) return 2;
        for (i = 2; i <= n; i += 2)
            result *= i;
    }else{
        if (n <= 1) return 1;
        for (i = 1; i <= n; i += 2)
            result *= i;
    }
    return result;
}

double fact_l_lambda(int l, int lambda)
{
    // 量子化学中册 P62 (10.6.5)
    return (double)factorial(l) / (factorial(lambda) * factorial(l - lambda));
}

double fi_l_ll_pax_pbx(int ii, int l1, int l2, double pax, double pbx, int flags)
{
// 《量子化学》中册 P63 第一个公式
    int i, j;
    double sum = 0;

    for (i = 0; i <= l1; i++) {
        for (j = 0; j <= l2; j++) {
            if (i + j == ii) {
                sum += (fact_l_lambda(l1, i) * fact_l_lambda(l2, j) * \
                                    pow(pax, l1 - i) * pow(pbx, l2 - j));
            }
        }
    }
    return  sum;
}

double I_xyz(int l1, double pax, int l2, double pbx, double gamma, int flags)
{
// 《量子化学》中册 P63 (10.6.8) (10.6.9) (10.6.10)
    int i;
    double sum = 0;

    for (i = 0; i <= (l1 + l2) / 2; i++) {
        if (flags == 1)
            printf("i = %d l1 = %d l2 = %d fi = %lf\n",i,l1,l2, fi_l_ll_pax_pbx(2*i, l1, l2, pax, pbx, flags));
        sum += factorial_2(2i - 1) / pow(2 * gamma, i) * sqrt(M_PI / gamma) * \
                fi_l_ll_pax_pbx(2*i, l1, l2, pax, pbx, flags);
    }
    return sum;
}

double gauss_K(double a, const gsl_vector *A, double b, const gsl_vector *B)
{
// 归一化系数 《量子化学》中册 P58 (10.4.1b)
// A, B 为坐标

    double result, norm_2;
    double x1, y1, z1, x2, y2, z2;

    x1 = gsl_vector_get(A, 0);
    y1 = gsl_vector_get(A, 1);
    z1 = gsl_vector_get(A, 2);

    x2 = gsl_vector_get(B, 0);
    y2 = gsl_vector_get(B, 1);
    z2 = gsl_vector_get(B, 2);


    norm_2 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2);
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
    double x1, x2, tmp;

    gsl_vector *center = gsl_vector_alloc(3);
    
    for (i = 0; i < 3; i++) {
        //x1 = gsl_vector_get(A, i);
        //x2 = gsl_vector_get(B, i);
        //tmp = (a * x1 +  b*x2) / gamma;
        //gsl_vector_set(center, i, tmp);
    // 与上面看似等同，但是计算的S与P轨道相互作用不一样
        center->data[i] = (a * A->data[i] + b * B->data[i]) / gamma;
    }

    // DEBUG
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

double overlap_gauss(const GTO g1, const gsl_vector* A, const GTO g2, const gsl_vector* B, int debug)
{
    double K, gamma, Ix, Iy, Iz, result = 0;
    double normal1, normal2;
    double coeff1, coeff2;
    gsl_vector *PA, *PB, *P;

    PA = gsl_vector_alloc(3);
    PB = gsl_vector_alloc(3);

    normal1 = g1.norm;
    normal2 = g2.norm;
    coeff1 = g1.coeff;
    coeff2 = g2.coeff;

    // compute the coordination of P. 《量子化学》中册, P77
    P = gaussian_product_center(g1.alpha, A, g2.alpha, B, debug);
    gsl_vector_memcpy(PA, A);
    gsl_vector_memcpy(PB, B);

    gsl_vector_sub(PA, P);
    gsl_vector_sub(PB, P);

    gamma = g1.alpha + g2.alpha;


    Ix = I_xyz(g1.l, -PA->data[0], g2.l, -PB->data[0], gamma, debug);
    Iy = I_xyz(g1.m, -PA->data[1], g2.m, -PB->data[1], gamma, debug);
    Iz = I_xyz(g1.n, -PA->data[2], g2.n, -PB->data[2], gamma, debug);

    K = gauss_K(g1.alpha, A, g2.alpha, B);

    //result += pow(M_PI/gamma, 1.5) * K * Ix * Iy * Iz * normal1 * normal2 * coeff1 * coeff2;
    result += K * Ix * Iy * Iz * normal1 * normal2 * coeff1 * coeff2;
    // doesn't do normalization
    // result += pow(M_PI/gamma, 1.5) * K * Ix * Iy * Iz;
    if (debug == 2) {
                printf("--------------------------------------------\n");
                vector_output(PA, 3, "PA:");
                vector_output(PB, 3, "PB:");
                gto_output(&g1, 1, "g1:");
                gto_output(&g2, 1, "g2:");
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
            result += overlap_gauss(g1, A, g2, B, debug);
        }
    }
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
