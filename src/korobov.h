/* C header for korobov.c *****************************************************/

#ifndef korobov_H
#define korobov_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

void korobov(int n, int d, int *generator, int randomize, double *res);
SEXP korobov_(SEXP n, SEXP d, SEXP generator, SEXP randomize);

#endif

