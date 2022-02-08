#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP cpwbart_it(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cpwbart_it", (DL_FUNC) &cpwbart_it, 3},
    {NULL, NULL, 0}
};

void R_init_BART(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

