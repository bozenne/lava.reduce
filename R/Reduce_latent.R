#' @title Specify the latent variables
#' @description Same as latent for latent variable model but check the consistency of the object
#' 
#' @param x \code{lvm}-object
#' @param value the names of the latent variables
`latent<-.lvm.reduced` <- function(x,value, ...){
 
  x <- lava:::`latent<-.lvm`(x, value = value, ...)
  
  if(!is.null(latent(x))){
    
    if(any(latent(x) %in% lp(x, type = "x"))){
      latent.pb <- latent(x)[latent(x) %in% lp(x, type = "x")]
      n.pb <- length(latent.pb)
      stop("the latent variable",if(n.pb>1){"s"}," \"",paste(latent.pb, collapse = "\" \""),"\" must not be contained in a linear predictor \n")
    }
  }
  
  return(x)
}
 