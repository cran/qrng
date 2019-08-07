/* C header for sobol.c *******************************************************/

#ifndef sobol_H
#define sobol_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define sobolMaxDim 16510
#define sobolMaxDegree 17
#define sobolMaxCol 32

void sobol(int n, int d, int randomize, double *res, int skip);
SEXP sobol_(SEXP n, SEXP d, SEXP randomize, SEXP skip);

#endif
