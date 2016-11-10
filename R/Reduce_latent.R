#' @title Specify the latent variables
#' @description Same as latent for latent variable model but check the consistency of the object
#' 
#' @param x \code{lvm}-object
#' @param value the names of the latent variables
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @export
`latent<-.lvm.reduced` <- function(x, ..., value){
 
  x <- callS3methodParent(x, FUN = "latent<-", class = "lvm.reduced", value = value, ...)
  
  if(!is.null(lava::latent(x))){
    
    if(any(lava::latent(x) %in% lp(x, type = "x"))){
      latent.pb <- lava::latent(x)[lava::latent(x) %in% lp(x, type = "x")]
      n.pb <- length(latent.pb)
      stop("the latent variable",if(n.pb>1){"s"}," \"",paste(latent.pb, collapse = "\" \""),"\" must not be contained in a linear predictor \n")
    }
  }
  
  return(x)
}
 