#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "hamiltonian.h"
#include "overlap.h"
#include "int2e.h"
#include "eri.h"
#include "int1e.h"

int main(int argc, char** argv)
{
    char *basis_base = NULL;
    int n;

    if (argc < 2)
        basis_base = "input_file";
    else
        basis_base = argv[1];

    INPUT_INFO *b = parse_input(basis_base);    
    n = b->basisCount;
    //printf("Basis Count: %d\n", n);

/*
    gsl_matrix* overlapMatrix = overlap_matrix(b);
    matrix_output(overlapMatrix, n, "OVERLAP INTEGRALS:");

    gsl_matrix* kmatrix = kinetic_matrix(b);
    matrix_output(kmatrix, n, "KINETIC INTEGRALS:");
    gsl_matrix* attractionmatrix = nuclear_attraction_matrix(b);
    matrix_output(attractionmatrix, n, "ATTRACTION INTEGRALS:");

    gsl_matrix* h = hamiltonian(b);
    matrix_output(h, n, "Hamiltonian Matirx:");
*/

#ifdef __INTEGRAL__INT2E__ONE__
    double *int2e;
#else
    double ****int2e;
#endif
    int2e = int2e_matrix(b);
    int2e_output(int2e, n, "TWO ELECTRON INTEGRAL:");


    return 0;
}
