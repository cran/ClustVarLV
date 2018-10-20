#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/
/* .Call calls */
extern SEXP _ClustVarLV_critcpp(SEXP, SEXP);
extern SEXP _ClustVarLV_mincpp(SEXP);
extern SEXP _ClustVarLV_powerEigen(SEXP);

static const R_CallMethodDef CallEntries[] = {
        {"_ClustVarLV_critcpp",    (DL_FUNC) &_ClustVarLV_critcpp,    2},
        {"_ClustVarLV_mincpp",     (DL_FUNC) &_ClustVarLV_mincpp,     1},
        {"_ClustVarLV_powerEigen", (DL_FUNC) &_ClustVarLV_powerEigen, 1},
        {NULL, NULL, 0}
};

void R_init_ClustVarLV(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);}

