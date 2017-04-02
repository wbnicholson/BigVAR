#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>



/* .Call calls */
extern SEXP BigVAR_ARFitVARXR(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_Eigencomp(SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_EigencompOO(SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_FistaElem(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_Fistapar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_gamloopElem(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_gamloopFista(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_GamLoopGL2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_GamLoopGLOO(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_gamloopHVAR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_gamloopOO(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_GamLoopSGL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_GamLoopSGLDP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_GamLoopSGLOO(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_GamLoopSGLOODP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_GamLoopSGLX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_GamLoopSGLXDP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_ICX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_norm2(SEXP);
extern SEXP BigVAR_powermethod(SEXP, SEXP);
extern SEXP BigVAR_proxvx2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_RelaxedLS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BigVAR_VARXCons(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"BigVAR_ARFitVARXR",     (DL_FUNC) &BigVAR_ARFitVARXR,      5},
    {"BigVAR_Eigencomp",      (DL_FUNC) &BigVAR_Eigencomp,       4},
    {"BigVAR_EigencompOO",    (DL_FUNC) &BigVAR_EigencompOO,     4},
    {"BigVAR_FistaElem",      (DL_FUNC) &BigVAR_FistaElem,       8},
    {"BigVAR_Fistapar",       (DL_FUNC) &BigVAR_Fistapar,        8},
    {"BigVAR_gamloopElem",    (DL_FUNC) &BigVAR_gamloopElem,    10},
    {"BigVAR_gamloopFista",   (DL_FUNC) &BigVAR_gamloopFista,   13},
    {"BigVAR_GamLoopGL2",     (DL_FUNC) &BigVAR_GamLoopGL2,     16},
    {"BigVAR_GamLoopGLOO",    (DL_FUNC) &BigVAR_GamLoopGLOO,    17},
    {"BigVAR_gamloopHVAR",    (DL_FUNC) &BigVAR_gamloopHVAR,    10},
    {"BigVAR_gamloopOO",      (DL_FUNC) &BigVAR_gamloopOO,      12},
    {"BigVAR_GamLoopSGL",     (DL_FUNC) &BigVAR_GamLoopSGL,     17},
    {"BigVAR_GamLoopSGLDP",   (DL_FUNC) &BigVAR_GamLoopSGLDP,   17},
    {"BigVAR_GamLoopSGLOO",   (DL_FUNC) &BigVAR_GamLoopSGLOO,   17},
    {"BigVAR_GamLoopSGLOODP", (DL_FUNC) &BigVAR_GamLoopSGLOODP, 17},
    {"BigVAR_GamLoopSGLX",    (DL_FUNC) &BigVAR_GamLoopSGLX,    17},
    {"BigVAR_GamLoopSGLXDP",  (DL_FUNC) &BigVAR_GamLoopSGLXDP,  17},
    {"BigVAR_ICX",            (DL_FUNC) &BigVAR_ICX,             8},
    {"BigVAR_norm2",          (DL_FUNC) &BigVAR_norm2,           1},
    {"BigVAR_powermethod",    (DL_FUNC) &BigVAR_powermethod,     2},
    {"BigVAR_proxvx2",        (DL_FUNC) &BigVAR_proxvx2,         6},
    {"BigVAR_RelaxedLS",      (DL_FUNC) &BigVAR_RelaxedLS,       6},
    {"BigVAR_VARXCons",       (DL_FUNC) &BigVAR_VARXCons,        8},
    {NULL, NULL, 0}
};

void R_init_BigVAR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

