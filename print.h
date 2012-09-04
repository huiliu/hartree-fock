#include <stdio.h>
#include "basis.h"
#include "common.h"

#ifndef INT_COMMON_PRINT__
#define INT_COMMON_PRINT__

// 以更好的格式输出矩阵
void matrix_output(const gsl_matrix *, int , char *);
// 以更好的格式输出向量
void vector_output(const gsl_vector *, int , char *);

void gto_output(const GTO* g, int count, char* msg);
void basis_set_output(const BASIS*, int, char* );
void atom_output(const ATOM_INFO** atom, int n);

#ifdef __INTEGRAL__INT2E__ONE__
    void int2e_output(double* e, int n, char* msg);
#else
    void int2e_output(double**** e, int n, char* msg);
#endif

#endif
