/* Register routines with R *****************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "korobov.h"
#include "ghalton.h"
#include "sobol.h"


static const R_CallMethodDef callMethods[] = {
	{"korobov_", (DL_FUNC) &korobov_, 4},
	{"ghalton_", (DL_FUNC) &ghalton_, 3},
	{"sobol_",   (DL_FUNC) &sobol_, 3},
	{NULL, NULL, 0}
};

void R_init_qrng(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL); /* s. WRE (2015, Section 5.4) */
}
