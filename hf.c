#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <assert.h>

//#define DEBUG_SCF
//#define DEBUG_s_root

void Load();
gsl_matrix* S_i_root(double *S, int n);
void matrix_output(gsl_matrix *m, int n, char *msg);
void vector_output(gsl_vector *v, int n, char *msg);
gsl_matrix *Fock(double h[][4], double e2_int[][4][4][4], double density[][4], int n);
void density(double Density[][4], gsl_matrix* coef, int n);
gsl_matrix* scf(gsl_matrix *f, const gsl_matrix *s_root, int n, double* energy);

int main(int argc, char** argv)
{
    Load();
    return 0;
}

void Load()
{
    int n = 4;
    int i, j, k, l;
    double e2[n][n][n][n];
    FILE *f;

    char *file_name = "2e";
    double S[] = {
                    0.100000000E+01, 0.658291970E+00, 0.460875162E+00, 0.511249073E+00,
                    0.658291970E+00, 0.100000000E+01, 0.511249073E+00, 0.856364449E+00,
                    0.460875162E+00, 0.511249073E+00, 0.100000000E+01, 0.658291970E+00,
                    0.511249073E+00, 0.856364449E+00, 0.658291970E+00, 0.100000000E+01
                   };
    double HH[][4] = {
    {-0.965997649E+00,-0.923907639E+00,-0.818582204E+00,-0.783996285E+00},
    {-0.923907639E+00,-0.928703434E+00,-0.783996285E+00,-0.857664629E+00},
    {-0.818582204E+00,-0.783996285E+00,-0.965997649E+00,-0.923907639E+00},
    {-0.783996285E+00,-0.857664629E+00,-0.923907639E+00,-0.928703434E+00}
                };
    double Density[4][4];

    char *file_coeff = "coeff";
    gsl_matrix *coeff  = gsl_matrix_alloc(n, n);

    //读入系数矩阵
    f = fopen(file_coeff, "r");
    gsl_matrix_fscanf(f, coeff);
    fclose(f);
    //matrix_output(coeff, n, "系数矩阵:\n");


    // 读入双电子积分
    f = fopen(file_name, "r");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++) {
                    //fscanf(f, "%lf", &e2[l][j][k][i]);
                    fscanf(f, "%lf", &e2[i][j][k][l]);
                }
            }
        }
    }
    fclose(f);
    i = 0;
    int itmax = 100;

    double energy, old_energy = 0.0;
    while (i < itmax) {
        printf("iter %d\n", i);
        //计算密度矩阵
        density(Density, coeff, n);
        // 计算FOCK矩阵
        gsl_matrix *F = Fock(HH, e2, Density, n);
        //matrix_output(F, n, "FOCK矩阵为:\n");
        gsl_matrix *s = S_i_root(S, n);
        coeff = scf(F, s, n, &energy);
        if (fabs(energy - old_energy) < 1.0E-6)
            break;
        else
            old_energy = energy;
        i++;
    }
}
void density(double Density[][4], gsl_matrix* coef, int n)
{
    int i, j, k;
    //printf("密度矩阵为:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Density[i][j] = 0;
            for (k = 0; k < 1; k++) {
                //Density[i][j] += 2 * coef[i][k] * coef[j][k];
                Density[i][j] += 2 * gsl_matrix_get(coef, i, k) * gsl_matrix_get(coef, j, k);
            }
            //printf("%20.10lf", Density[i][j]);
        }
        //printf("\n");
    }
}

gsl_matrix* scf(gsl_matrix *f, const gsl_matrix *s_root, int n, double* energy)
{
    gsl_matrix *eigVector = gsl_matrix_alloc(n, n);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(2*n);
    gsl_vector *eigValue = gsl_vector_alloc(n);

    gsl_matrix *ft = gsl_matrix_calloc(n, n);
    gsl_matrix *ftt = gsl_matrix_calloc(n, n);
    
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, s_root, f, 1.0, ft);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ft, s_root, 1.0, ftt);

    //求本征矢量和本征值
    //gsl_eigen_symm(b, dialg_S, w);
    gsl_eigen_symmv(ftt, eigValue, eigVector, w);
    *energy = gsl_vector_min(eigValue);

    vector_output(eigValue, n, "FOCK本征值为：\n");
#ifdef DEBUG_SCF
    matrix_output(eigVector, n,  "FOCK本征矢量为:\n");
#endif

    gsl_matrix *c = gsl_matrix_calloc(n, n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, s_root, eigVector, 1.0, c);
#ifdef DEBUG_SCF
    matrix_output(c, n,  "新的轨道系数为:\n");
#endif
    return c;
}

// 计算Fock矩阵
gsl_matrix *Fock(double h[][4], double e2_int[][4][4][4], double density[][4], int n)
{
    int i, j, k, l;
    double tmp = 0;
    gsl_matrix *m = gsl_matrix_alloc(n,n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            tmp = 0;
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++)
                    tmp += ((2*e2_int[i][j][k][l] - e2_int[i][l][k][j]) * density[l][k]/2);
            }
            gsl_matrix_set(m, i, j, tmp + h[i][j]);
        }
    }
    return m;
}

// 此函数用于计算重叠矩阵逆阵的平方根
/*
 *  P * A * P^(-1) = D
 *
 *  sqrt(A) = P^(-1) * sqrt(D) * P
 *
 */
gsl_matrix* S_i_root(double *S, int n)
{

    gsl_matrix_view gv = gsl_matrix_view_array(S, n, n);
    gsl_matrix *a = gsl_matrix_alloc_from_matrix(&gv.matrix,0, 0, n, n);
    gsl_matrix *b = gsl_matrix_alloc(n, n);
    gsl_matrix *p = gsl_matrix_alloc(n, n);
    
    gsl_matrix_memcpy(b, a);
    gsl_matrix_free(a);
#ifdef DEBUG_s_root
    matrix_output(b, n, "初始重叠矩阵为:\n");
#endif
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(2*n);
    gsl_vector *dialg_S = gsl_vector_alloc(n);

    //求本征矢量和本征值
    //gsl_eigen_symm(b, dialg_S, w);
    gsl_eigen_symmv(b, dialg_S, p, w);
#ifdef DEBUG_s_root
    vector_output(dialg_S, n, "重叠矩阵本征值为：\n");
    matrix_output(p, n,  "重叠矩阵本征矢量为:\n");
#endif

    // 将本征值开方
    int i;
    gsl_matrix *s_root = gsl_matrix_calloc(n, n);
    for (i = 0; i < n; i++) {
        gsl_matrix_set(s_root, i, i, 1/sqrt(dialg_S->data[i]));
        //gsl_matrix_set(s_root, i, i, dialg_S->data[i]);
    }
    //matrix_output(s_root, n, "the root of S^-1:\n");


    // 利用LU分解求本征矢量的逆
    gsl_matrix *pp = gsl_matrix_alloc(n, n);
    gsl_matrix *inverse = gsl_matrix_alloc(n, n);
    gsl_permutation *permutation = gsl_permutation_alloc(n);
    int s;

    gsl_matrix_memcpy(pp, p);

    gsl_linalg_LU_decomp(pp, permutation, &s);
    gsl_linalg_LU_invert(pp, permutation, inverse);
#ifdef DEBUG_s_root
    matrix_output(inverse, n, "本征矢量的逆:\n");
#endif
    gsl_matrix_free(pp);
    gsl_permutation_free(permutation);

    // 两个矩阵相乘，测试逆阵计算是否正确
    //gsl_matrix *c = gsl_matrix_calloc(n, n);
    //gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, p, inverse, 1.0, c);
    //matrix_output(c, n, "test:\n");
    //gsl_matrix_free(c);
    
    // 求得逆阵的平方根
    gsl_matrix *c = gsl_matrix_calloc(n, n);
    gsl_matrix *d = gsl_matrix_calloc(n, n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, p, s_root, 1.0, c);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, c, inverse, 1.0, d);
#ifdef DEBUG_s_root
    matrix_output(d, n, "逆阵的平方根:\n");
#endif

    gsl_matrix_free(b);
    gsl_matrix_free(p);
    //gsl_matrix_free(d);
    return d;
}

void matrix_output(gsl_matrix *m, int n, char *msg)
{
    int i, j;

    printf("%s",msg);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%20.11g", gsl_matrix_get(m, i, j));
        printf("\n");
    }
}

void vector_output(gsl_vector *v, int n, char *msg)
{
    int i;

    printf("%s",msg);
    for (i = 0; i < n; i++)
        printf("%15.08g", gsl_vector_get(v, i));
    printf("\n");
}
