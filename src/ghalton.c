/* C function for computing a generalized Halton sequence *********************/

#include "ghalton.h"


/**
 * @title Generate n Points of a d-dimensional Generalized Halton Sequence
 * @param n number of points
 * @param d dimension
 * @param method character string indicating which sequence is generated
 *        (generalized Halton or (plain) Halton)
 * @param res pointer to the result matrix
 * @return void
 * @author Marius Hofert based on C. Lemieux's RandQMC
 */
void ghalton(int n, int d, const char *method, double *res)
{
        static int perm[ghaltonMaxDim];
        int base, i, j, k, l, count, maxindex, f;
        double u, U;
        unsigned int tmp;
        unsigned int shcoeff[ghaltonMaxDim][32]; /* the coefficients of the shift */
        unsigned int coeff[32];

        /* Init */
        GetRNGstate();
        for(j=0; j<d; j++) {
                base = primes[j];
                u = 0;
                for(k=31; k >= 0; k--) {
                        U = unif_rand();
                        shcoeff[j][k] = (int) (base * U);
                        u += shcoeff[j][k];
                        u /= base;
                }
                res[j*n] = u;
        }
        PutRNGstate();

        /* Main */
	Rboolean is_ghalton;
        if(strcmp(method, "generalized") == 0) is_ghalton=TRUE; else is_ghalton=FALSE;
        if(!is_ghalton) {
                for(j=0;j<d;j++) { perm[j] = 1; }
        } else {
                for(j=0; j<d; j++) { perm[j] = permTN2[j]; }
        }
        for(i=1; i<n; i++) {
                for(j=0; j<d; j++) {
                        tmp = i;
                        base = primes[j]; /* (j+1)st prime number for this dimension */
                        memset (&coeff, 0, sizeof(int) * 32); /* clear the coefficients */

                        /* Find i in the prime base */
                        k = 0;
                        while((tmp > 0) && (k < 32)) {
                               coeff[k] = tmp % base;
                               tmp /= base;
                               k++;
                        }
                        maxindex = k;
                        for(l=maxindex + 1; l<32; l++){ coeff[l] = 0; }
                        u = 0.0;
                        k = 31;
                        f = perm[j];
                        while(k >= 0) {
                                u += (f * coeff[k] + shcoeff[j][k]) % base;
                                u /= base;
                                k--;
                        }
                        res[j*n + i] = u;
                }
        }
}

/**
 * @title R Interface to C for Generating a Generalized Halton Sequence
 * @param n number of points
 * @param d dimension
 * @param method character string indicating which sequence is generated
 *        (generalized Halton or (plain) Halton)
 * @return (n, d)-matrix
 * @author Marius Hofert
 */
SEXP ghalton_(SEXP n, SEXP d, SEXP method)
{
    /* Input parameters */
    int n_ = asInteger(n); /* numeric(1) */
    int d_ = asInteger(d); /* numeric(1) */
    const char *method_ = CHAR(STRING_ELT(method, 0)); /* character(1) */

    /* Create result object */
    SEXP res = PROTECT(allocMatrix(REALSXP, n_, d_)); /* (n,d)-matrix */
    double *res_ = REAL(res); /* pointer to the values of res */

    /* Main */
    ghalton(n_, d_, method_, res_);

    /* Return */
    UNPROTECT(1); /* clean-up */
    return res;
}
