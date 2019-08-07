/* C header for ghalton.c *****************************************************/

#ifndef ghalton_H
#define ghalton_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define ghaltonMaxDim 360

void ghalton(int n, int d, const char *method, double *res);
SEXP ghalton_(SEXP n, SEXP d, SEXP method);

#endif
