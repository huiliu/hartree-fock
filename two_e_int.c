#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "overlap.h"
#include "basis.h"

/*
int B(int i1, int i2, int r1, int r2, double pax, double pbx, double qcx, double qdx)
{
    double result;
    for (i = 0; i <
        for (i = 0; i <
            for (i = 0; i <)
                tmp = pow(-1, i2) * fi_l_ll_pax_pbx(l1, l2, pax, pbx) * fi_l_ll_pax_pbx(l3, l4, qcx, qdx);
                tmp *= (factorial(i1) * factorial(i2) / (pow(4*gamma1, i1) * pow(4*gamma2, i2) * pow(delta, i1+i2)))
                tmp *= ((pow(4*gamma1, r1) * pow(4*gamma2, r2) * pow(delta, r1+r2)) / (factorial(r1) * factorial(r2) * \
                        factorial(i1 - 2*r1) * factorial(i2 - 2*r2)));
                tmp *= ((factorial(i1 + i2 - 2*(r1+r2)) * pow(-1, u) * pow(px, i1 + i2 - 2*(r1+r2) - 2*u) * \
                        power(delta, u)) / (factorial(u) * factorial(i1 + i2 - 2*(r1+r2-u))))
}
*/
void main(void)
{
    int i, j, k, l, ii, jj;
    double JK[basis_count][basis_count][basis_count][basis_count];

    for (i = 0; i < basis_count; i++) {
        for (j = 0; j < basis_count; j++) {
            for (k = 0; k < basis_count; k++) {
                for (l = 0; l < basis_count; l++) {

                    for (ii = 0; ii < gauss_count; ii++) {
                        for (jj = 0; jj < gauss_count; jj++) {
                            sum
}

double coulomb_repulsion_trans(const gsl_vector* A, const GTO* g1, \
                         const gsl_vector* B, const GTO* g2, \
                         const gsl_vector* C, const GTO* g3, \
                         const gsl_vector* D, const GTO* g4)
{
    double xa, xb, xc, xd, ya, yb, yc, yd, za, zb, zc, zd;
    double alphaa, alphab, alphac, alphad;
    double la, lb, lc, ld, ma, mb, mc, md, na, nb, nc, nd;
    double sum = 0;

    sum = coulomb_repulsion(xa, ya, za, norma, la, ma, na, alphaa, \
                            xb, yb, zb, normb, lb, mb, nb, alphab, \
                            xc, yc, zc, normc, lc, mc, nc, alphac, \
                            xd, yd, zd, normd, ld, md, nd, alphad)

    return sum;

    /*
    double rab2, rcd2,rpq2,xp,yp,zp,xq,yq,zq,gamma1,gamma2,delta,sum;

    gamma1 = g1->alpha + g2->alpha;
    gamma2 = g3->alpha + g4->alpha;

    gsl_vector *AB = gsl_vector_alloc(3);
    gsl_vector *CD = gsl_vector_alloc(3);
    gsl_vector *PQ = gsl_vector_alloc(3);
    gsl_vector_memcpy(AB, A);
    gsl_vector_memcpy(CD, C);
    double rab2 = pow(gsl_blas_dnrm2(AB), 2);
    double rcd2 = pow(gsl_blas_dnrm2(CD), 2);

    gsl_vector* P = gaussian_product_center(g1->alpha, A, g2->alpha, B);
    gsl_vector* Q = gaussian_product_center(g3->alpha, C, g4->alpha, D);
    gsl_vector_memcpy(PQ, P);
    double rpq2 = pow(gsl_blas_dnrm2(PQ), 2);
    
    delta = 0.25 * (1/gamma1 + 1/gamma2);
    */

}
