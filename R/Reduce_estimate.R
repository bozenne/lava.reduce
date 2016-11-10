#' @title Estimate a latent variable model with linear predictors (LP)
#' @name estimate.lvm.reduced
#' @description Add columns corresponding to the LPs in the dataset (filled with 0) so that LPs are not treated as a latent variable
#' 
#' @param x \code{lvm.reduced}-object
#' @param data \code{data.frame}
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @example 
#' R/examples/EX_reducedModel.R
#' @export
estimate.lvm.reduced <- function(x, data, ...){
  
  ####  add columns
  name.LP <- lp(x)
  data <- cbind(data,
                matrix(0,nrow = NROW(data), ncol = length(name.LP), dimnames = list(NULL,name.LP))
  )
  
  return(callS3methodParent(x, FUN = "estimate", class = "lvm.reduced", data=data, ...))
}