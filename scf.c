#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "common.h"
#include "scf.h"
#include "int2e.h"
#include "basis.h"
#include "overlap.h"
#include "hamiltonian.h"

//#define DEBUG_SCF
//#define DEBUG_s_root

int main(int argc, char** argv)
{
    char *inputFile, *coeffFile;
    if (argc < 2) {
        inputFile = "input_file";
        coeffFile = NULL;
    }else if (argc == 2) {
        inputFile = argv[1];
        coeffFile = NULL;
    }else if (argc == 3) {
        inputFile = argv[1];
        coeffFile = argv[2];
    }

    HartreeFock(inputFile, coeffFile);
    return 0;
}

void HartreeFock(char* fname, char* coeffFile)
{
    int i, n;
    FILE *f;
    double** Density;
    INPUT_INFO *b;

    b = parse_input(fname);

    n = b->basisCount;
    //printf("TOTAL ELECTRON: %d\n", b->eCount);
    gsl_matrix* S = overlap_matrix(b);
    gsl_matrix* H = hamiltonian(b);
#ifdef __INTEGRAL__INT2E__ONE__
    double *int2e;
#else
    double ****int2e;
#endif
    int2e = int2e_matrix(b);

// ---------- READ COEFF DATA FROM FILE --------
    Density = (double**)malloc(sizeof(double*)*n);
    for (i = 0; i < n; i++)
        *(Density+i) = malloc(sizeof(double)*n);

    char *file_coeff = "coeff";
    gsl_matrix *coeff  = gsl_matrix_alloc(n, n);

    if (coeffFile == NULL)
        file_coeff = "coeff";
    else
        file_coeff = coeffFile;
    //读入系数矩阵
    f = fopen(file_coeff, "r");
    gsl_matrix_fscanf(f, coeff);
    fclose(f);
// -----------------END --------------------

    i = 0;
    int itmax = 200;

    double energy, old_energy = 0.0;

    for (i = 0; i < itmax; i++) {
        printf("\niter %d\n", i);
        //计算密度矩阵
        density(Density, coeff, b->eCount, n);
        // 计算FOCK矩阵
        gsl_matrix *F = Fock(H, int2e, Density, n);
#ifdef DEBUG_SCF
        matrix_output(F, n, "FOCK矩阵为:");
#endif

        gsl_matrix *s = S_i_root(S, n);
        coeff = scf(F, s, n, &energy);
        printf("E = %15.8lf%15.8lf\n", old_energy, energy);
        if (fabs(energy - old_energy) < 1.0E-6)
            break;
        else
            old_energy = energy;
    }
}
void density(double **Density, gsl_matrix* coef, int e, int n)
{
    int i, j, k;

#ifdef DEBUG_SCF
    printf("%s %d\n","DENSITY MATRIX:", e);
#endif
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Density[i][j] = 0;
            for (k = 0; k < e/2; k++) {
                Density[i][j] += 2 * gsl_matrix_get(coef, i, k) * gsl_matrix_get(coef, j, k);
            }
#ifdef DEBUG_SCF
            printf("%15.6E", Density[i][j]);
#endif
        }
#ifdef DEBUG_SCF
        printf("\n");
#endif
    }
#ifdef DEBUG_SCF
    gsl_matrix_free(coef);
#endif
}

gsl_matrix* scf(gsl_matrix *f, const gsl_matrix *s_root, int n, double* energy)
{
    gsl_matrix *c = gsl_matrix_calloc(n, n);
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
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eigValue, eigVector, GSL_EIGEN_SORT_VAL_ASC);
    //gsl_eigen_symmv_sort(eigValue, eigVector, GSL_EIGEN_SORT_VAL_DESC);
    *energy = gsl_vector_min(eigValue);

    vector_output(eigValue, n, "FOCK本征值为:");
#ifdef DEBUG_SCF
    matrix_output(eigVector, n,  "FOCK本征矢量为:");
    matrix_output(s_root, n, "SQRT of Overlap:");
#endif

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, s_root, eigVector, 1.0, c);
#ifdef DEBUG_SCF
    matrix_output(c, n,  "新的轨道系数为:");
#endif
    return c;
}

/* 构造Fock矩阵
 *  h       hamiltonian matrix(one electron integrals).
 *  e2_int  two electron integrals.
 *  density density matrix.
 *  n       count of basis function.
 */
gsl_matrix *Fock(gsl_matrix *h, double ****e2_int, double** density, int n)
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
            gsl_matrix_set(m, i, j, tmp + gsl_matrix_get(h, i, j));
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
gsl_matrix* S_i_root(gsl_matrix *S, int n)
{

    gsl_matrix *a = gsl_matrix_alloc_from_matrix(S,0, 0, n, n);
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
    gsl_eigen_symmv(b, dialg_S, p, w);
    gsl_eigen_symmv_free(w);
#ifdef DEBUG_s_root
    vector_output(dialg_S, n, "重叠矩阵本征值为：\n");
    matrix_output(p, n,  "重叠矩阵本征矢量为:\n");
#endif

    // 将本征值开方
    int i;
    gsl_matrix *s_root = gsl_matrix_calloc(n, n);
    for (i = 0; i < n; i++) {
        gsl_matrix_set(s_root, i, i, 1/sqrt(dialg_S->data[i]));
    }


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
