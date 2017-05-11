#' @docType package
#' @title lava.reduce
#' @name lava.reduce
#' @description Latent variable models with linear predictors. Should enable a faster optimisation when dealing with a large number of variables. 
#' @import lava
#' @importFrom mets loglikMVN scoreMVN
#' @importFrom numDeriv jacobian
#' @importFrom stats as.formula coef delete.response predict setNames terms
#' @importFrom utils getS3method
#' @useDynLib lava.reduce, .registration=TRUE
NULL


#####' @importFrom lava addhook endogenous estimate exogenous index intercept latent manifest ordinal regression vars "latent<-"
