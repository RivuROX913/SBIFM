// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// updateEta_c
arma::mat updateEta_c(arma::mat Lambda, arma::vec ps, int k, arma::mat Y, int n);
RcppExport SEXP _SBIFM_updateEta_c(SEXP LambdaSEXP, SEXP psSEXP, SEXP kSEXP, SEXP YSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ps(psSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(updateEta_c(Lambda, ps, k, Y, n));
    return rcpp_result_gen;
END_RCPP
}
// updateLambda_c
arma::mat updateLambda_c(arma::mat eta, arma::mat Plam, arma::vec ps, arma::mat Y, int k, int p);
RcppExport SEXP _SBIFM_updateLambda_c(SEXP etaSEXP, SEXP PlamSEXP, SEXP psSEXP, SEXP YSEXP, SEXP kSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Plam(PlamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ps(psSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(updateLambda_c(eta, Plam, ps, Y, k, p));
    return rcpp_result_gen;
END_RCPP
}
// updatePsi_c
arma::mat updatePsi_c(double df, arma::mat Lambda, arma::vec tauh, int p, int k);
RcppExport SEXP _SBIFM_updatePsi_c(SEXP dfSEXP, SEXP LambdaSEXP, SEXP tauhSEXP, SEXP pSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tauh(tauhSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(updatePsi_c(df, Lambda, tauh, p, k));
    return rcpp_result_gen;
END_RCPP
}
// updateDeltaTauh_c
List updateDeltaTauh_c(arma::mat Lambda, arma::mat psijh, double ad1, int p, int k, double bd1, arma::vec delta, arma::vec tauh, double ad2, double bd2);
RcppExport SEXP _SBIFM_updateDeltaTauh_c(SEXP LambdaSEXP, SEXP psijhSEXP, SEXP ad1SEXP, SEXP pSEXP, SEXP kSEXP, SEXP bd1SEXP, SEXP deltaSEXP, SEXP tauhSEXP, SEXP ad2SEXP, SEXP bd2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type psijh(psijhSEXP);
    Rcpp::traits::input_parameter< double >::type ad1(ad1SEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type bd1(bd1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tauh(tauhSEXP);
    Rcpp::traits::input_parameter< double >::type ad2(ad2SEXP);
    Rcpp::traits::input_parameter< double >::type bd2(bd2SEXP);
    rcpp_result_gen = Rcpp::wrap(updateDeltaTauh_c(Lambda, psijh, ad1, p, k, bd1, delta, tauh, ad2, bd2));
    return rcpp_result_gen;
END_RCPP
}
// updateSigma_c
arma::vec updateSigma_c(arma::vec tmp2, int p, double as, int n, double bs);
RcppExport SEXP _SBIFM_updateSigma_c(SEXP tmp2SEXP, SEXP pSEXP, SEXP asSEXP, SEXP nSEXP, SEXP bsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type tmp2(tmp2SEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    rcpp_result_gen = Rcpp::wrap(updateSigma_c(tmp2, p, as, n, bs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SBIFM_updateEta_c", (DL_FUNC) &_SBIFM_updateEta_c, 5},
    {"_SBIFM_updateLambda_c", (DL_FUNC) &_SBIFM_updateLambda_c, 6},
    {"_SBIFM_updatePsi_c", (DL_FUNC) &_SBIFM_updatePsi_c, 5},
    {"_SBIFM_updateDeltaTauh_c", (DL_FUNC) &_SBIFM_updateDeltaTauh_c, 10},
    {"_SBIFM_updateSigma_c", (DL_FUNC) &_SBIFM_updateSigma_c, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_SBIFM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}