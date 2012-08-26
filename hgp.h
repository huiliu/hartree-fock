#include "overlap.h"
#include "common.h"
#include "basis.h"


#ifndef INTEGRALS_HGP__ERI
#define INTEGRALS_HGP__ERI

#define GetShellBasisCount(i, L)   \
    if (i == 0) L = 1; /* S */\
    else if (i == 1) L = 3; /* P */\
    else if (i == 2) L = 6; /* D */\
    else if (i == 3) L = 10; /* F */\
    else if (i == 4) L = 15; /* G */\
    else if (i == 5) L = 21; /* H */\
    else{ \
        fprintf(stderr, "This Program cann't deal with the orbital of L >5\n"); \
        exit(EXIT_FAILURE); \
    }
        
void HGPShell(double ****, const BASIS *, 
                        int *, int *, int *, int *, int, int );
double HGPBasisHRR(BASIS *, BASIS *, BASIS *b3, BASIS *,
                 gsl_vector *, gsl_vector *, double *, int );
double HGPBasis(const BASIS* b1, const BASIS* b2,
                const BASIS* b3, const BASIS* b4,
                const gsl_vector *AB, const gsl_vector *CD,
                double *XSXS, int debug);
double HGPHrrVRR(int l1, int m1, int n1, int l3, int m3, int n3,
            double zeta, double gamma, double ro,
            const gsl_vector *PA, const gsl_vector *PB, const gsl_vector *QC,
            const gsl_vector *QD, const gsl_vector *WP, const gsl_vector *WQ,
            int m, double *T);
int HGPIndex(int, int, int, int, int, int);
/*
#define HGPIndex(l1, m1, n1, N1, l3, m3, n3, N2)    gsl_pow_5(MAXSHELL, l1) + \
                gsl_pow_4(MAXSHELL, m1) + gsl_pow_3(MAXSHELL, n1) + \
                gsl_pow_2(MAXSHELL, l3) + MAXSHELL * m3 + n3 
*/
#endif
