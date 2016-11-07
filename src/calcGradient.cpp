//// NOTE: http://thread.gmane.org/gmane.comp.lang.r.rcpp/4457

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <Rmath.h>
// #include "mets.h" // example biglasso calls bigmemory
// #include "mets/mvtdst.h"


using namespace Rcpp ;
using namespace std ;
using namespace arma ;


void calcLP(arma::mat& data, const arma::vec& p, unsigned n_lp,
            const std::vector<IntegerVector >& indexCoef, const std::vector<IntegerVector >& indexEndo, 
            const IntegerVector& indexLP);

// [[Rcpp::export]]
arma::mat scoreLVM(arma::mat data, const arma::vec& p, 
                   Function scoreFun, const arma::mat& mu, const arma::mat& dmu, const arma::mat& S, const arma::mat& dS,
                   const std::vector<IntegerVector >& indexCoef, const std::vector<IntegerVector >& indexEndo,
                   const IntegerVector& indexIntercept,
                   const IntegerVector& indexLP, const IntegerVector& indexManifest,
                   bool indiv){
  
  unsigned n_lp = indexLP.size(), n_endo;
  
  //// linear predictor
  calcLP(data, p, n_lp, indexCoef, indexEndo, indexLP);
  
  //// score
  arma::mat Y = data.cols(as<uvec>(indexManifest));
  arma::mat score = as<mat>(scoreFun(Y, mu, dmu, S, dS));
  
  // SEXP A = loglikMVN(Y, // yl
  //                    NULL, // yu
  //                    NULL, // status
  //                    mu,dmu, S, dS, // 
  //                    NULL, NULL, NULL, NULL, NULL, true); // z su dsu threshold dthreshold score
                     
  
  //// chain rule
  arma::mat X;
  
  for(unsigned iterLP=0 ; iterLP < n_lp ; iterLP++){
    n_endo = indexEndo[iterLP].size();
    for(unsigned iterEndo=0 ; iterEndo < n_endo ; iterEndo++){
      score.col(indexCoef[iterLP][iterEndo]) = data.col(indexEndo[iterLP][iterEndo]) % score.col(indexIntercept[iterLP]);
      // score = X %*% score(intercept) where %*% is the Schur product, i.e. element wise
    }
  }
  
  if(indiv == false){
    score = sum(score,0);
  }
  
  return(score);
}

// other functions

void calcLP(arma::mat& data, const arma::vec& p, unsigned n_lp,
            const std::vector<IntegerVector >& indexCoef, const std::vector<IntegerVector >& indexEndo, const IntegerVector& indexLP){
  
  arma::colvec coef, LP;
  arma::mat X;

  for(unsigned iterLP=0 ; iterLP < n_lp ; iterLP++){

    coef = p.elem(as<uvec>(indexCoef[iterLP]));
    X = data.cols(as<uvec>(indexEndo[iterLP]));
    data.col(indexLP[iterLP]) = X * coef;

  }
  
}

