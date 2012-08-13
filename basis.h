#include <gsl/gsl_vector.h>

#ifndef __INTEGRAL__BASIS__
#define __INTEGRAL__BASIS__

// angstrom to bohr radius
#define ANGS_2_BOHR     1.88971616463
#define ATOM_MAX_COUNT  10

typedef struct gto {
    short int l, m, n;
    double alpha;
    double coeff;
    double norm;
}GTO;

typedef struct _b {
// the parameter of basis function
    short int gaussCount;
    GTO* gaussian;   // 3 表示一个基函数由3个gaussian函数构成
    gsl_vector* xyz;
}BASIS;

typedef struct atom_INFORMATON_ {
// 定义一个原子信息, 包括核电荷数，基函数数目，元素符号，坐标，基函数
    short int n;              // Atomic Number
    short int basisCount;    // basis function count
    char symbol[3];     // Element Symbol
    gsl_vector* coordination;
}ATOM_INFO;

typedef struct FILE_INPUT_ {
// 存储解析后的整个输入文件
    short int atomCount;          // atom count in current system
    short int eCount;             // electron count
    short int basisCount;    // basis function count
    BASIS* basisSet;   // save all of the basis set
    gsl_vector *gXYZ[ATOM_MAX_COUNT];   // save all of coordination
    ATOM_INFO *atomList[ATOM_MAX_COUNT];   // element information
}INPUT_INFO;

void* read_basis(const char * );

INPUT_INFO* parse_input(const char* file_name);
BASIS* readbasis(FILE * f, int basisCount);
double normalize_coeff(const GTO *);
void bridge(BASIS* b, ATOM_INFO** atomList, int atomCount);
int gtoIsNeg(const GTO* g);
#endif
