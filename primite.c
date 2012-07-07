/*
 * 此文件的目的主要是为了存储基函数，
 *
 * 当前以《量子化学》徐光宪，中册，P240计算NH3为例子
 * 存储其中的基函数，并实现乘法
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG2

// 声明结构用于保存Gaussian函数的参数
typedef struct g{
    double coef[4];
}GAUSS;

typedef struct gto {
    double a;
    double A;
}GTO;

// 因为使用的NH3一个基函数用三个Gaussian函数拟合
typedef struct N2_basis{
    GAUSS gauss[3];
    struct N2_basis *next;
}BASIS;

typedef struct Int {
    double coeff;
    double expon;
}INT_PARAM;

BASIS* Primitive(char *fname)
{
// 循环读入基函数
// 问题：
//      怎么合适的返回数值
//
    FILE *f;
    int n, i, j;
    BASIS *b, *bhead;

    // HEAD为一个空数据
    bhead = b = malloc(sizeof(BASIS));
    f = fopen(fname, "r");
    while (fscanf(f, "%d", &n) != EOF) {
        b->next = malloc(sizeof(BASIS));
        b = b->next;
        // n 表示每行的Gaussian函数数目
        for (i = 0; i < n; i++) {
            // 一个Gaussian函数由四个参数构成
            for (j = 0; j < 4; j++)
                fscanf(f, "%lf", &b->gauss[i].coef[j]);
        }
        b->next = NULL;
    }
    //b->next = malloc(sizeof(BASIS));
    //b->next->next=NULL;
    fclose(f);
    return bhead;
}

// 打印基函数参数，n 为基函数数目
void output(const BASIS* p)
{
    BASIS *ptr = (BASIS *)p;
    //BASIS *ptr = p->next;
    int i, j;
    while ((ptr = ptr->next) != NULL) {
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 4; j++)
                printf("%20.10lf", ptr->gauss[i].coef[j]);
            printf("\n");
        }
        printf("\n");   //打印完一个基函数
    }
}

BASIS* multi(const BASIS *p_basis)
{
    int i, j;
    int ii = 0, jj = 0;
    double coeff, expon, Int_Value;
    BASIS *p_basis_a = (BASIS*)p_basis;
    BASIS *p_basis_b = (BASIS*)p_basis;

    double Int_Result[2][2];

    while ((p_basis_a = p_basis_a->next) != NULL) {
        p_basis_b = (BASIS*)p_basis;
        jj = 0;
        while ((p_basis_b = p_basis_b->next) != NULL) {
            if (ii == jj) {
                Int_Value = 0;
                for (i = 0; i < 3; i++) {
                    for (j = 0; j < 3; j++) {
                        // 计算指数前面的系数
                        coeff = p_basis_a->gauss[i].coef[0]*p_basis_a->gauss[i].coef[1]* \
                            p_basis_b->gauss[j].coef[0]*p_basis_b->gauss[j].coef[1];
                        expon = p_basis_a->gauss[i].coef[2] + p_basis_b->gauss[j].coef[2];
                        printf("%20.10lf%20.10lf\n", coeff, expon);
                        Int_Value += (4 * M_PI * coeff * 0.25 * sqrt(M_PI / pow(expon, 3)));
                        //printf("%20.10lf\n", Int_Value);
                    }
                    printf("\n");
                }
                printf("积分值为:\t%20.10lf\n\n", Int_Value);
                Int_Result[ii][jj] = Int_Value;
            }
            jj++;
        }
        ii++;
    }
    
    return NULL;
}

int factorial(int n)
{
    int i, result = 1;

    if (n <= 1)
        return 1;

    for (i = 1; i <= n; i++)
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

    for (i = 0; i < l; i++) {
        for (j = 0; j < ll; j++) {
            if ((i + j) == ii)
                sum += (fact_l_lambda(l, i) * fact_l_lambda(ll, j) * pow(pax, l - i) * pow(pbx, ll - j));        
        }
    }

    return  sum;
}

double I_x(double a, double b, int l, int ll, double pax, double pbx)
{
// 《量子化学》中册 P63 (10.6.8) (10.6.9) (10.6.10)
    int i;
    double sum = 0;

    for (i = 0; i < (l + ll) / 2; i++)
        sum += fi_l_ll_pax_pbx(i, l, ll, pax, pbx) * factorial(factorial(2i - 1)) / pow(2 * (a + b), i);
    return sum;
}

double normal_coeff(double a, double b, const gsl_vector *A, const gsl_vector *B)
{
// 归一化系数 《量子化学》中册 P58 (10.4.1b)
    // A, B 为坐标
    gsl_vector *AB = gsl_vector_alloc(3);
    gsl_vector_memcpy(AB, A);
    gsl_vector_sub(AB, B);

    double result = exp(-a * b / (a + b) * pow(gsl_blas_dnrm2(AB), 2));
    gsl_vector_free(AB);

    return result;
}

int main(int argc, char **argv)
{
    BASIS *basis;

    basis = Primitive("p");
    output(basis);
    multi(basis);
    return 0;
}
