#' @title Hook for \code{lvm.reduced} objects
#' @name hookLVMreduced
#' @description Pre and post processing methods
#' 
#' @inheritParams estimate.lvm.reduced
#' @param estimator the function used to compute the derivatives of the objective function
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @export
lavaReduce.estimate.hook <- function(x,data,estimator,...) {

    dots <- list(...)
  if("lvm.reduced" %in% class(x) && length(lp(x))>0){

    ## set reference endogeneous variable where there is no linear predictor
    if(length(latent(x))>0){
      n.latent <- length(latent(x))
      for(iterLatent in 1:n.latent){
        endoLatent <- names(which(x$index$M[lava::latent(x)[iterLatent],]==1))
        reference <- setdiff(endoLatent, names(x$lp))
        if(length(reference)==0){
          stop("A reduce LVM must have for each latent variable at least one observed outcome not in reduced form \n")
        }else{
          f <- stats::as.formula(paste(reference, lava::latent(x)[iterLatent], sep = "~"))
          lava::regression(x, f) <- 1
          f <- stats::as.formula(paste(reference, 1, sep = "~"))
          lava::intercept(x, f) <- 0
        }
      }
    }
    
      ## 
      if(is.null(dots$optim$start)){ # intialisation of the parameters
          ## non LP
          x0 <- kill(x, value = lp(x), restaure = TRUE) # remove LP and external parameters
          x0 <- clean(x0)
		
          startLVM <- initializer.lavaReduce(x0, data = data, optim = dots$optim)
        
        ## update
        dots$optim$start <- stats::setNames(rep(0, length = length(coef(x))), coef(x))
        dots$optim$start[names(startLVM)] <- startLVM
        dots$optim$start <- dots$optim$start[which(!is.na(dots$optim$start))]
    }
    
    ## update estimator
    if(estimator == "gaussian"){
      estimator <- paste0("gaussian", lava.options()$estimator.default.reduce)
    }
   
    validEstimator <- paste0("gaussian",c("",1:3))
    if(estimator %in% validEstimator){
      estimator <- paste0(estimator,"LP")
    } else {
      stop("reduced estimator for ",estimator," not implemented \n",
           "available estimators: ",paste(validEstimator, collapse = " "),"\n")
    }
    
    
  }       
    return( c(list(x=x, data=data, estimator = estimator),dots) )
}

#' @rdname hookLVMreduced
#' @export
lavaReduce.post.hook <- function(x){
  
  return(x)
  
}
