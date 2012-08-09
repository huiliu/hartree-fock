#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "common.h"
#include "hamiltonian.h"
#include "overlap.h"
#include "int2e.h"
#include "eri.h"
#include "int1e.h"

int main(int argc, char** argv)
{
    int n, m;
    char *basis_base = NULL, randchar[3];

    if (argc < 2) {
        fprintf(stderr, "%s\n", "Please give a input file.");
        exit(EXIT_FAILURE);
    }else
        basis_base = argv[1];

    m = strlen(basis_base);
    FMW *fmw;
#ifdef ERI_INT_USE_REDBLACK_TREE
    int i;
    CALLOC(fmw, sizeof(FMW), 1);
    CALLOC(fmw->RBList, sizeof(struct rbtree), 36);
    fmw->count = 36;
    for (i = 0; i < fmw->count; i++)
        fmw->RBList[i] = *rbinit(compare, NULL);
#else
    fmw = NULL;
#endif

    INPUT_INFO *b = parse_input(basis_base);    
    n = b->basisCount;
    //printf("Basis Count: %d\n", n);

/*
    gsl_matrix* overlapMatrix = overlap_matrix(b);
    matrix_output(overlapMatrix, n, basis_base, "ov-");

    gsl_matrix* kmatrix = kinetic_matrix(b);
    matrix_output(kmatrix, n, basis_base, "km-");

    gsl_matrix* attractionmatrix = nuclear_attraction_matrix(b, fmw);
    matrix_output(attractionmatrix, n, basis_base, "at-");

    gsl_matrix* h = hamiltonian(b, fmw);
    matrix_output(h, n, basis_base, "hm-");
*/

#ifdef __INTEGRAL__INT2E__ONE__
    double *int2e;
#else
    double ****int2e;
#endif
    int2e = int2e_matrix(b, fmw);
    //int2e_output(int2e, n, basis_base, "2e-");

    return 0;
}
