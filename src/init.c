/* GENERATED BY THE R FUNCTION CALL:
 * tools::package_native_routine_registration_skeleton("gcdnet", character_only = FALSE)
 */

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(loofit)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(traverse_zigzag)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sunicoldstart)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(suniwalpha)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lslassonet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(soft_unilassonet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"loofit",        (DL_FUNC) &F77_NAME(loofit),        8},
    {"traverse_zigzag", (DL_FUNC) &F77_NAME(traverse_zigzag), 7},
    {"sunicoldstart", (DL_FUNC) &F77_NAME(sunicoldstart), 28},
    {"suniwalpha",    (DL_FUNC) &F77_NAME(suniwalpha),    28},
    {"lslassonet",    (DL_FUNC) &F77_NAME(lslassonet),    25},
    {"soft_unilassonet", (DL_FUNC) &F77_NAME(soft_unilassonet), 26},
    {NULL, NULL, 0}
};

void R_init_gcdnet(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
