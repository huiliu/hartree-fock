#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "basis.h"
#include "eri_drive.h"
#include "eri_os.h"

double**** ERI_Matrix(INPUT_INFO* b)
{
/*
    compute eri matrix element
 */
    int debug = 0;
    register int i, j, k, l;
    int basis_count, atomCount;
    int is_dup = 0;
    BASIS *basisSet;

    //ATOM_INFO **alist = b->atomList;
    basis_count = b->basisCount;
    basisSet = b->basisSet;
    atomCount = b->atomCount;

    // use four dimension array output int2e
    double ****e2;
    e2 = (double****)Malloc(sizeof(double***)*basis_count);
    for (i = 0; i < basis_count; i++) {
        *(e2+i) = (double***)Malloc(sizeof(double**)*basis_count);
        for (j = 0; j < basis_count; j++) {
            *(e2[i]+j) = (double**)Malloc(sizeof(double*)*basis_count);
            for (k = 0; k < basis_count; k++) {
                *(e2[i][j]+k) = (double*)Calloc(sizeof(double), basis_count);
            }
        }
    }

    //omp_set_num_threads(2);
    //#pragma omp parallel for private(j, k, l)
    for (i = 0; i < basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            for (k = 0; k < basis_count; k++) {
                for (l = 0; l < basis_count; l++) {
                    //debug = 0;
                    //if (i == 0 && j == 2 && k == 0 && l == 2)
                    //if (i == 2 && j == 0 && k == 0 && l == 2)
                    //    debug = 999;

                    ChkERISym(e2, i, j, k, l, basis_count, &is_dup);
                    if (is_dup)     continue;

//            fprintf(stdout, "---%d---%d---%d---%d---\n", i, j, k, l);
                    e2[i][j][k][l] = \
                    e2[i][j][l][k] = \
                    e2[j][i][k][l] = \
                    e2[j][i][l][k] = \
                    e2[k][l][i][j] = \
                    e2[k][l][j][i] = \
                    e2[l][k][i][j] = \
                    e2[l][k][j][i] = ERI_basis_OS(&basisSet[i], 
                                                 &basisSet[j], 
                                                 &basisSet[k], 
                                                 &basisSet[l], 
                                                 NULL, debug);
                }
            }
        }
    }
    return e2;
}
