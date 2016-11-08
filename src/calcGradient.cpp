// [[Rcpp::depends(RcppArmadillo, mets)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <Rmath.h>
#include <mets.h> 

using namespace Rcpp ;
using namespace std ;
using namespace arma ;


void calcLP(arma::mat& data, const arma::vec& p, unsigned n_lp,
            const std::vector<IntegerVector >& indexCoef, const std::vector<IntegerVector >& indexEndo, 
            const IntegerVector& indexLP);

// [[Rcpp::export]]
arma::mat scoreLVM(arma::mat data, const arma::vec& p, 
                   const arma::mat& mu, const SEXP& dmu, const arma::mat& S, const SEXP& dS,
                   const std::vector<IntegerVector >& indexCoef, const std::vector<IntegerVector >& indexEndo,
                   const IntegerVector& indexIntercept,
                   const IntegerVector& indexLP, const IntegerVector& indexManifest,
                   bool indiv){
  
  unsigned n_lp = indexLP.size(), n_endo;
  
  //// linear predictor
  calcLP(data, p, n_lp, indexCoef, indexEndo, indexLP);
  
  //// score
  arma::mat Y = data.cols(as<uvec>(indexManifest));
  arma::mat score = mets::_loglikMVN(Y, // yl
                                     R_NilValue, // yu
                                     R_NilValue, // status
                                     mu, dmu, S, dS, //
                                     R_NilValue, R_NilValue, R_NilValue, R_NilValue, R_NilValue, true); // z su dsu threshold dthreshold score
  
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



