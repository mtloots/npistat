#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern void C_npi_bern(const int*, const int*, const int*, const int*, const int*, double*, double*);
extern void C_npi_mfuture(const int*, const double*, const double*, const int*, const int*,
                          const int*, const int*, const double*, double*, double*);
extern void C_npi_untie(const int*, const double*, const double*, const int*, const double*,
                        const int*, const double*, double*, double*);
extern void C_npi_share(const int*, const double*, const double*, const int*, double*);
extern void C_npi_place(const int*, const double*, const double*, const int*, const double*,
                        double*, double*, double*);
extern void C_npi_mbounds(const int*, const double*, const double*, const double*, const int*,
                          const int*, const double*, double*, double*);
extern void C_npi_turnbull(const int*, const double*, const double*, const int*, const double*,
                           int*, double*, double*, double*);
extern void C_npi_bounds(const int*, const double*, const double*, const int*, const double*,
                         const int*, const int*, const double*, double*, double*, double*, double*,
                         double*, double*, double*, double*, double*, double*);
extern void C_npi_rc(const int*, const double*, const int*, double*, double*, double*);

static const R_CMethodDef CEntries[] = {
    {"C_npi_bern",     (DL_FUNC) &C_npi_bern,      7},
    {"C_npi_mfuture",  (DL_FUNC) &C_npi_mfuture,  10},
    {"C_npi_untie",    (DL_FUNC) &C_npi_untie,     9},
    {"C_npi_share",    (DL_FUNC) &C_npi_share,     5},
    {"C_npi_place",    (DL_FUNC) &C_npi_place,     8},
    {"C_npi_mbounds",  (DL_FUNC) &C_npi_mbounds,   9},
    {"C_npi_turnbull", (DL_FUNC) &C_npi_turnbull,  9},
    {"C_npi_bounds",   (DL_FUNC) &C_npi_bounds,   18},
    {"C_npi_rc",       (DL_FUNC) &C_npi_rc,        6},
    {NULL, NULL, 0}
};

void R_init_npistat(DllInfo *dll){
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
