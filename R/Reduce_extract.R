# {{{ doc
#' @title Extract variables from a reduced latent variable model
#' @name getReduce
#' @description Extract variables as in a standard lvm but including or not the variables that compose the linear predictor
#' 
#' @param x \code{lvm}-object
#' @param lp should the name of the linear predictors be returned?
#' @param xlp should the name of the variables that the linear predictors aggregates be returned?
#' @param top for compatibility with endogenous.lvm
#' @param latent logical defining whether latent variables without parents should be included in the result 
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @details lp returns all the linear predictors of the \code{lvm}-object. 
#' The other functions plays the same role as those defined in the lava package.
#' 
#' @examples  
#' ## regression
#' m <- lvm.reduced()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' vars(m)
#' vars(m, lp = TRUE)
#' vars(m, lp = FALSE, xlp = TRUE)
#' 
#' exogenous(m)
#' exogenous(m, lp = FALSE)
#' exogenous(m, xlp = TRUE)
#' 
#' endogenous(m)
#' endogenous(m, lp = FALSE) # should not change
#' endogenous(m, xlp = TRUE) # should not change
#' 
#' coef(m)
#' 
#' ## lvm
#' m <- lvm.reduced()
#' m <- regression(m, x=paste0("x",1:10),y="y1", reduce = TRUE)
#' m <- regression(m, x=paste0("x",51:150),y="y2", reduce = TRUE)
#' covariance(m) <- y1~y2
#' 
#' vars(m)
#' vars(m, lp = FALSE)
#' 
#' exogenous(m)
#' exogenous(m, lp = FALSE)
#' 
#' endogenous(m)
#' endogenous(m, lp = FALSE) # should not change
# }}}


# {{{ vars.lvm.reduced
#' @rdname getReduce 
#' @export
vars.lvm.reduced <- function(x, lp = TRUE, xlp = FALSE, ...){
  if(xlp){
    hiddenX <- lp(x, type = "x")
  }else{
    hiddenX <- NULL
  }
  if(lp){
    names.lp <- NULL
  }else{
    names.lp <- lp(x, type = "name", ...)
  }
  
  class(x) <- setdiff(class(x),"lvm.reduced")
  allVars <- unique(c(vars(x), hiddenX))
  
  return(setdiff(allVars,names.lp))
}
# }}}

# {{{ vars.lvmfit.reduced
#' @rdname getReduce 
#' @export
vars.lvmfit.reduced <- function(x, lp = FALSE, xlp = TRUE, ...){
  return(vars(Model(x), lp = lp, xlp = xlp, ...))
}
# }}}

# {{{ exogenous.lvm.reduced
#' @rdname getReduce 
#' @export
exogenous.lvm.reduced <- function(x, lp = TRUE, xlp = FALSE, ...){
  
  if(xlp){
    hiddenX <- lp(x, type = "x")
  }else{
    hiddenX <- NULL
  }
  if(lp){
    names.lp <- NULL
  }else{
    names.lp <- lp(x, type = "name")
  }
  
  class(x) <- setdiff(class(x),"lvm.reduced")
  allExo <- unique(c(exogenous(x), hiddenX))
  
  return(setdiff(allExo, names.lp))
}
# }}}

# {{{ exogenous.lvmfit.reduced
#' @rdname getReduce 
#' @export
exogenous.lvmfit.reduced <- function(x, lp = FALSE, xlp = TRUE, ...){
  return(exogenous(Model(x), lp = lp, xlp = xlp, ...))
}
# }}}

# {{{ endogenous.lvm.reduced
#' @rdname getReduce 
#' @export
endogenous.lvm.reduced <- function(x, top = FALSE, latent = FALSE, ...){
  
  observed <- lava::manifest(x, lp = FALSE, xlp = FALSE)
  if (latent) observed <- vars(x, lp = FALSE, xlp = FALSE)
  if (top) {
    M <- x$M
    res <- c()
    for (i in observed)
      if (!any(M[i,]==1))
        res <- c(res, i)
    return(res)
  }
  exo <- exogenous(x, lp = FALSE, xlp = FALSE)
  return(setdiff(observed,exo))
  
}
# }}}

# {{{ endogenous.lvmfit.reduced
#' @rdname getReduce 
#' @export
endogenous.lvmfit.reduced <- function(x, lp = FALSE, xlp = TRUE, ...){
  return(endogenous(Model(x), lp = lp, xlp = xlp, ...))
}
# }}}

# {{{ manifest.lvm.reduced
#' @rdname getReduce
#' @export
manifest.lvm.reduced <- function(x, lp = TRUE, xlp = FALSE, ...) {
  
  vars <- vars(x, lp = lp, xlp = xlp)
  if (length(vars)>0)
    setdiff(vars, lava::latent(x))
  else
    NULL
}
# }}}

# {{{ manifest.lvmfit.reduced
#' @rdname getReduce 
#' @export
manifest.lvmfit.reduced <- function(x, lp = FALSE, xlp = TRUE, ...){
  return(manifest(Model(x), lp = lp, xlp = xlp, ...))
}
# }}}
