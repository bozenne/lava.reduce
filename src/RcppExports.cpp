// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// scoreLVM
arma::mat scoreLVM(arma::mat data, const arma::vec& p, const arma::mat& mu, const SEXP& dmu, const arma::mat& S, const SEXP& dS, const std::vector<IntegerVector >& indexCoef, const std::vector<IntegerVector >& indexEndo, const IntegerVector& indexIntercept, const IntegerVector& indexLP, const IntegerVector& indexManifest, bool indiv);
RcppExport SEXP lava_reduce_scoreLVM(SEXP dataSEXP, SEXP pSEXP, SEXP muSEXP, SEXP dmuSEXP, SEXP SSEXP, SEXP dSSEXP, SEXP indexCoefSEXP, SEXP indexEndoSEXP, SEXP indexInterceptSEXP, SEXP indexLPSEXP, SEXP indexManifestSEXP, SEXP indivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type dmu(dmuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type dS(dSSEXP);
    Rcpp::traits::input_parameter< const std::vector<IntegerVector >& >::type indexCoef(indexCoefSEXP);
    Rcpp::traits::input_parameter< const std::vector<IntegerVector >& >::type indexEndo(indexEndoSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type indexIntercept(indexInterceptSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type indexLP(indexLPSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type indexManifest(indexManifestSEXP);
    Rcpp::traits::input_parameter< bool >::type indiv(indivSEXP);
    rcpp_result_gen = Rcpp::wrap(scoreLVM(data, p, mu, dmu, S, dS, indexCoef, indexEndo, indexIntercept, indexLP, indexManifest, indiv));
    return rcpp_result_gen;
END_RCPP
}
