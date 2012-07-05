/*
 * 此文件的目的主要是为了存储基函数，
 *
 * 当前以《量子化学》徐光宪，中册，P240计算NH3为例子
 * 存储其中的基函数，并实现乘法
 *
 */
#include <stdio.h>
#include <stdlib.h>

#define DEBUG2

// 声明结构用于保存Gaussian函数的参数
typedef struct g{
    double coef[4];
}GAUSS;

// 因为使用的NH3一个基函数用三个Gaussian函数拟合
typedef struct N2_basis{
    GAUSS gauss[3];
    struct N2_basis *next;
}BASIS;

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
    double coeff, expon;
    BASIS *p_basis_a = (BASIS*)p_basis;
    BASIS *p_basis_b = (BASIS*)p_basis;

    while ((p_basis_a = p_basis_a->next) != NULL) {
        p_basis_b = (BASIS*)p_basis;
        while ((p_basis_b = p_basis_b->next) != NULL) {
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    // 计算指数前面的系数
                    coeff = p_basis_a->gauss[i].coef[0]*p_basis_a->gauss[i].coef[1]* \
                        p_basis_b->gauss[j].coef[0]*p_basis_b->gauss[j].coef[1];
                    expon = p_basis_a->gauss[i].coef[2] + p_basis_b->gauss[j].coef[2];
                    printf("%20.10lf%20.10lf\n", coeff, expon);
                }
                printf("\n");
            }
        }
    }
    
    return NULL;
}

int main(int argc, char **argv)
{
    BASIS *basis;

    basis = Primitive("p");
    output(basis);
    multi(basis);
    return 0;
}
