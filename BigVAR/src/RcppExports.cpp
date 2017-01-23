// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/BigVAR.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// VARXCons
MatrixXd VARXCons(NumericMatrix Y1, NumericMatrix X1, const int k, const int p, const int m, int s, bool oos, bool contemp);
RcppExport SEXP BigVAR_VARXCons(SEXP Y1SEXP, SEXP X1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP mSEXP, SEXP sSEXP, SEXP oosSEXP, SEXP contempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< bool >::type oos(oosSEXP);
    Rcpp::traits::input_parameter< bool >::type contemp(contempSEXP);
    rcpp_result_gen = Rcpp::wrap(VARXCons(Y1, X1, k, p, m, s, oos, contemp));
    return rcpp_result_gen;
END_RCPP
}
// ARFitVARXR
List ARFitVARXR(NumericMatrix K21, const int k, const int p, int m, int s);
RcppExport SEXP BigVAR_ARFitVARXR(SEXP K21SEXP, SEXP kSEXP, SEXP pSEXP, SEXP mSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type K21(K21SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(ARFitVARXR(K21, k, p, m, s));
    return rcpp_result_gen;
END_RCPP
}
// ICX
List ICX(NumericMatrix Y1, NumericMatrix X1, double k, int pmax, int smax, double m, std::string pen, int h);
RcppExport SEXP BigVAR_ICX(SEXP Y1SEXP, SEXP X1SEXP, SEXP kSEXP, SEXP pmaxSEXP, SEXP smaxSEXP, SEXP mSEXP, SEXP penSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pmax(pmaxSEXP);
    Rcpp::traits::input_parameter< int >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< std::string >::type pen(penSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(ICX(Y1, X1, k, pmax, smax, m, pen, h));
    return rcpp_result_gen;
END_RCPP
}
// gamloopFista
cube gamloopFista(NumericVector beta_, const mat& Y, const mat& Z, const colvec gammgrid, const double eps, const colvec& YMean2, const colvec& ZMean2, mat& B1, int k, int p, double tk, int k1, int s);
RcppExport SEXP BigVAR_gamloopFista(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP tkSEXP, SEXP k1SEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const colvec >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type tk(tkSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopFista(beta_, Y, Z, gammgrid, eps, YMean2, ZMean2, B1, k, p, tk, k1, s));
    return rcpp_result_gen;
END_RCPP
}
// Eigencomp
List Eigencomp(mat& Z1, List groups, int n1, int k1);
RcppExport SEXP BigVAR_Eigencomp(SEXP Z1SEXP, SEXP groupsSEXP, SEXP n1SEXP, SEXP k1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    rcpp_result_gen = Rcpp::wrap(Eigencomp(Z1, groups, n1, k1));
    return rcpp_result_gen;
END_RCPP
}
// EigencompOO
List EigencompOO(mat& ZZ1, List groups, int n1, int k);
RcppExport SEXP BigVAR_EigencompOO(SEXP ZZ1SEXP, SEXP groupsSEXP, SEXP n1SEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type ZZ1(ZZ1SEXP);
    Rcpp::traits::input_parameter< List >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(EigencompOO(ZZ1, groups, n1, k));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopGL2
List GamLoopGL2(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y1, const mat& Z1, List jj, List jjfull, List jjcomp, double eps, const colvec& YMean2, const colvec& ZMean2, int k, int pk, const List M2f_, const List eigvalF_, const List eigvecF_);
RcppExport SEXP BigVAR_GamLoopGL2(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigvalF_SEXP, SEXP eigvecF_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const List >::type eigvalF_(eigvalF_SEXP);
    Rcpp::traits::input_parameter< const List >::type eigvecF_(eigvecF_SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopGL2(beta_, Activeset, gamm, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M2f_, eigvalF_, eigvecF_));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopGLOO
List GamLoopGLOO(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y, const mat& Z, List jj, List jjfull, List jjcomp, double eps, colvec& YMean2, colvec& ZMean2, int k, int pk, List M2f_, List eigvalF_, List eigvecF_, int k1);
RcppExport SEXP BigVAR_GamLoopGLOO(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigvalF_SEXP, SEXP eigvecF_SEXP, SEXP k1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< List >::type eigvalF_(eigvalF_SEXP);
    Rcpp::traits::input_parameter< List >::type eigvecF_(eigvecF_SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopGLOO(beta_, Activeset, gamm, Y, Z, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M2f_, eigvalF_, eigvecF_, k1));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLOO
List GamLoopSGLOO(NumericVector beta_, const List Activeset_, const NumericVector gamm, const double alpha, const mat& Y, const mat& Z, List jj_, const List jjfull_, List jjcomp_, const double eps, const colvec& YMean2, const colvec& ZMean2, const int k1, const int pk, const List M2f_, const NumericVector eigs_, double m);
RcppExport SEXP BigVAR_GamLoopSGLOO(SEXP beta_SEXP, SEXP Activeset_SEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP jj_SEXP, SEXP jjfull_SEXP, SEXP jjcomp_SEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP k1SEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigs_SEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const List >::type Activeset_(Activeset_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< List >::type jj_(jj_SEXP);
    Rcpp::traits::input_parameter< const List >::type jjfull_(jjfull_SEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp_(jjcomp_SEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< const int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< const int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type eigs_(eigs_SEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLOO(beta_, Activeset_, gamm, alpha, Y, Z, jj_, jjfull_, jjcomp_, eps, YMean2, ZMean2, k1, pk, M2f_, eigs_, m));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLOODP
List GamLoopSGLOODP(NumericVector beta_, const List Activeset_, mat gamm, const colvec alpha, const mat& Y, const mat& Z, List jj_, const List jjfull_, List jjcomp_, const double eps, const colvec& YMean2, const colvec& ZMean2, const int k1, const int pk, const List M2f_, const NumericVector eigs_, double m);
RcppExport SEXP BigVAR_GamLoopSGLOODP(SEXP beta_SEXP, SEXP Activeset_SEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP jj_SEXP, SEXP jjfull_SEXP, SEXP jjcomp_SEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP k1SEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigs_SEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const List >::type Activeset_(Activeset_SEXP);
    Rcpp::traits::input_parameter< mat >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< List >::type jj_(jj_SEXP);
    Rcpp::traits::input_parameter< const List >::type jjfull_(jjfull_SEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp_(jjcomp_SEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< const int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< const int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type eigs_(eigs_SEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLOODP(beta_, Activeset_, gamm, alpha, Y, Z, jj_, jjfull_, jjcomp_, eps, YMean2, ZMean2, k1, pk, M2f_, eigs_, m));
    return rcpp_result_gen;
END_RCPP
}
// Fistapar
mat Fistapar(const mat Y, const mat Z, const mat phi, const int L, const double lambda, const double eps, const double tk, const int k);
RcppExport SEXP BigVAR_Fistapar(SEXP YSEXP, SEXP ZSEXP, SEXP phiSEXP, SEXP LSEXP, SEXP lambdaSEXP, SEXP epsSEXP, SEXP tkSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const mat >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double >::type tk(tkSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(Fistapar(Y, Z, phi, L, lambda, eps, tk, k));
    return rcpp_result_gen;
END_RCPP
}
// gamloopHVAR
cube gamloopHVAR(NumericVector beta_, const mat& Y, const mat& Z, colvec gammgrid, const double eps, const colvec& YMean2, const colvec& ZMean2, mat& B1, const int k, const int p);
RcppExport SEXP BigVAR_gamloopHVAR(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< colvec >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopHVAR(beta_, Y, Z, gammgrid, eps, YMean2, ZMean2, B1, k, p));
    return rcpp_result_gen;
END_RCPP
}
// gamloopOO
cube gamloopOO(NumericVector beta_, const mat Y, const mat Z, colvec gammgrid, const double eps, const colvec YMean2, const colvec ZMean2, mat B1, const int k, const int p, colvec w, List groups_);
RcppExport SEXP BigVAR_gamloopOO(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP wSEXP, SEXP groups_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< colvec >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< colvec >::type w(wSEXP);
    Rcpp::traits::input_parameter< List >::type groups_(groups_SEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopOO(beta_, Y, Z, gammgrid, eps, YMean2, ZMean2, B1, k, p, w, groups_));
    return rcpp_result_gen;
END_RCPP
}
// FistaElem
mat FistaElem(const mat& Y, const mat& Z, mat phi, const int p, const int k, double lambda, const double eps, const double tk);
RcppExport SEXP BigVAR_FistaElem(SEXP YSEXP, SEXP ZSEXP, SEXP phiSEXP, SEXP pSEXP, SEXP kSEXP, SEXP lambdaSEXP, SEXP epsSEXP, SEXP tkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< mat >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double >::type tk(tkSEXP);
    rcpp_result_gen = Rcpp::wrap(FistaElem(Y, Z, phi, p, k, lambda, eps, tk));
    return rcpp_result_gen;
END_RCPP
}
// gamloopElem
cube gamloopElem(NumericVector beta_, const mat& Y, const mat& Z, colvec gammgrid, const double eps, const colvec YMean2, const colvec ZMean2, mat B1, const int k, const int p);
RcppExport SEXP BigVAR_gamloopElem(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< colvec >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopElem(beta_, Y, Z, gammgrid, eps, YMean2, ZMean2, B1, k, p));
    return rcpp_result_gen;
END_RCPP
}
// powermethod
List powermethod(mat A, colvec x1);
RcppExport SEXP BigVAR_powermethod(SEXP ASEXP, SEXP x1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< colvec >::type x1(x1SEXP);
    rcpp_result_gen = Rcpp::wrap(powermethod(A, x1));
    return rcpp_result_gen;
END_RCPP
}
// norm2
double norm2(NumericVector x);
RcppExport SEXP BigVAR_norm2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(norm2(x));
    return rcpp_result_gen;
END_RCPP
}
// RelaxedLS
mat RelaxedLS(const mat K, mat B2, int k, int p, int k1, int s);
RcppExport SEXP BigVAR_RelaxedLS(SEXP KSEXP, SEXP B2SEXP, SEXP kSEXP, SEXP pSEXP, SEXP k1SEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< mat >::type B2(B2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(RelaxedLS(K, B2, k, p, k1, s));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLX
List GamLoopSGLX(NumericVector beta_, List Activeset, NumericVector gamm, double alpha, const mat& Y1, const mat& Z1, List jj, List jjfull, List jjcomp, double eps, colvec YMean2, colvec ZMean2, int k, int pk, List M2f_, NumericVector eigs, int k1);
RcppExport SEXP BigVAR_GamLoopSGLX(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigsSEXP, SEXP k1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eigs(eigsSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLX(beta_, Activeset, gamm, alpha, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M2f_, eigs, k1));
    return rcpp_result_gen;
END_RCPP
}
// proxvx2
colvec proxvx2(colvec v2, int L, double lambda, int m, int k, int F1);
RcppExport SEXP BigVAR_proxvx2(SEXP v2SEXP, SEXP LSEXP, SEXP lambdaSEXP, SEXP mSEXP, SEXP kSEXP, SEXP F1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< colvec >::type v2(v2SEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type F1(F1SEXP);
    rcpp_result_gen = Rcpp::wrap(proxvx2(v2, L, lambda, m, k, F1));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGL
List GamLoopSGL(NumericVector beta_, List Activeset, const NumericVector gamm, const double alpha, const mat& Y1, const mat& Z1, List jj, const List jjfull, const List jjcomp, const double eps, const colvec YMean2, const colvec ZMean2, const int k, const int pk, const List M1f_, const List M2f_, const NumericVector eigs_);
RcppExport SEXP BigVAR_GamLoopSGL(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M1f_SEXP, SEXP M2f_SEXP, SEXP eigs_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< const List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< const List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M1f_(M1f_SEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type eigs_(eigs_SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGL(beta_, Activeset, gamm, alpha, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M1f_, M2f_, eigs_));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLDP
List GamLoopSGLDP(NumericVector beta_, List Activeset, const mat gamm, const colvec alpha, const mat& Y1, const mat& Z1, List jj, const List jjfull, const List jjcomp, const double eps, const colvec YMean2, const colvec ZMean2, const int k, const int pk, const List M1f_, const List M2f_, const NumericVector eigs_);
RcppExport SEXP BigVAR_GamLoopSGLDP(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M1f_SEXP, SEXP M2f_SEXP, SEXP eigs_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< const mat >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< const List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< const List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M1f_(M1f_SEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type eigs_(eigs_SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLDP(beta_, Activeset, gamm, alpha, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M1f_, M2f_, eigs_));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLXDP
List GamLoopSGLXDP(NumericVector beta_, List Activeset, mat gamm, colvec alpha, const mat& Y1, const mat& Z1, List jj, List jjfull, List jjcomp, double eps, colvec YMean2, colvec ZMean2, int k, int pk, List M2f_, NumericVector eigs, int k1);
RcppExport SEXP BigVAR_GamLoopSGLXDP(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigsSEXP, SEXP k1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< mat >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eigs(eigsSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLXDP(beta_, Activeset, gamm, alpha, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M2f_, eigs, k1));
    return rcpp_result_gen;
END_RCPP
}
