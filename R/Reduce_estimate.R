#' @title Estimate a latent variable model with linear predictors (LP)
#' @description Add columns corresponding to the LPs in the dataset (filled with 0) so that LPs are not treated as a latent variable
#' 
#' @example 
#' R/examples/EX_reducedModel.R
#' @export
estimate.lvm.reduced <- function(x, data, ...){
  
  ####  add columns
  name.LP <- lp(x)
  data <- cbind(data,
                data.frame(matrix(0,nrow = NROW(data), ncol = length(name.LP), dimnames = list(NULL,name.LP)))
  )
  return(estimate.lvm(x, data=data, ...))
  #return(lava:::estimate.lvm(x, data=data, ...))
}