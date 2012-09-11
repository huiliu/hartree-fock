#include <gsl/gsl_vector.h>

#ifndef __INTEGRAL__BASIS__
#define __INTEGRAL__BASIS__

// angstrom to bohr radius
#define ANGS_2_BOHR     1.88971616463

typedef struct {
    double x, y, z;
}COORD;

typedef struct gto {
    unsigned short gtoID;
    short int l, m, n;
    double alpha;
    double coeff;
}GTO;

typedef struct _b {
// the parameter of basis function
    unsigned short int gaussCount;
    unsigned short int L;
    int l, m, n;
    GTO* gaussian;   // 3 表示一个基函数由3个gaussian函数构成
    gsl_vector *xyz;
}BASIS;

typedef struct atom_INFORMATON_ {
// 定义一个原子信息, 包括核电荷数，基函数数目，元素符号，坐标，基函数
    unsigned short int n;                   // Atomic Number
    unsigned short int basisCount;          // basis function count
    char symbol[3];                         // Element Symbol
    gsl_vector *coordination;
}ATOM_INFO;

typedef struct FILE_INPUT_ {
// 存储解析后的整个输入文件
    unsigned short int atomCount;           // atom count in current system
    unsigned short int eCount;              // electron count
    unsigned short int basisCount;          // basis function count
    unsigned short int gtoCount;            // basis function count
    BASIS *basisSet;                        // save all of the basis set
    GTO *gtoSet;
    COORD **P;
    gsl_vector **gXYZ;                      // save all of coordination
    ATOM_INFO** atomList;                   // element information
}INPUT_INFO;

INPUT_INFO* parse_input(const char* file_name);
BASIS* readbasis(FILE *f, int, unsigned short *);
double normalize_coeff(const GTO *);
void bridge(BASIS* b, ATOM_INFO** atomList, int atomCount);
int gtoIsNeg(const GTO* g);
COORD **PList(const BASIS *b, int n, int gtoCount);



#endif
