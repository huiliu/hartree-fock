#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#ifndef __INTEGRAL__COMMON__
#define __INTEGRAL__COMMON__
// 以更好的格式输出矩阵
void matrix_output(const gsl_matrix *, int , char *);
// 以更好的格式输出向量
void vector_output(const gsl_vector *, int , char *);
#endif
