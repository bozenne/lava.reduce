#' @docType package
#' @title lavaReduce
#' @name lavaReduce
#' @description Latent variable models with linear predictors. Should enable a faster optimisation when dealing with a large number of variables. 
#' @import lava
#' @importFrom mets loglikMVN scoreMVN
#' @importFrom numDeriv jacobian
#' @importFrom stats as.formula coef delete.response model.frame predict setNames terms
#' @importFrom utils getS3method
#' @useDynLib lavaReduce, .registration=TRUE
NULL


#####' @importFrom lava addhook endogenous estimate exogenous index intercept latent manifest ordinal regression vars "latent<-"
