# {{{ lavaReduce.estimate.hook
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
# }}}

# {{{ lavaReduce.post.hook
#' @rdname hookLVMreduced
#' @export
lavaReduce.post.hook <- function(x){
  
  return(x)
  
}
# }}}

# {{{ lavaReduce.clean.hook
#' @rdname clean
#' @export
lavaReduce.clean.hook <- function(x, rm.exo,
                                   rm.lp = TRUE, simplify.reduce = TRUE, simplify, ...){

    if("lvm.reduced" %in% class(x)){
        if(!missing(simplify)){
            simplify.reduce <- simplify
        }

        ## remove empty lp
        n.link <- lp(x, type = "n.link")
        if(rm.lp && any(n.link==0)){ 
            name.lp <- names(n.link)[which(n.link==0)]
            x <- kill(x, value = name.lp)
            x$lp[lp(x, lp = name.lp, type = "endo")] <- NULL
        }

        ## remove useless exogeneous variables
        if(length(lp(x, type = "link", format = "vector")) == 0 ){
            ## lvm.reduced object with no linear predictor
            if(simplify.reduce){
                class(x) <- setdiff(class(x), "lvm.reduced")
            }
        }else{
            ## lvm.reduced object with linear predictors
            var.exogenous <- exogenous(x)
            M.reg <- x$index$A
            M.cov <- x$index$P - diag(diag(x$index$P))

            if(rm.exo && length(var.exogenous) > 0){
                indexClean.reg <- which(rowSums(M.reg[var.exogenous,,drop = FALSE]!=0)==0)
                indexClean.cov <- which(rowSums(M.cov[var.exogenous,,drop = FALSE]!=0)==0)
                #indexClean.lp <- which(var.exogenous %in% lp(x, type = "x") == FALSE)
                #indexClean <- intersect(indexClean.reg, intersect(indexClean.cov, indexClean.lp))
                indexClean <- intersect(indexClean.reg, indexClean.cov)
                if(length(indexClean)>0){
                    x <- kill(x, value =  var.exogenous[indexClean])
                }
            }

            rm.exo <- FALSE
         
        }
    }

    return(list(x=x,
                rm.exo=rm.exo))

}
# }}}

# {{{ Description kill and cancel 
#' @title Remove variables from linear predictors
#' @description Remove one or several variable from the linear predictors
#' @name kill.lp
#' 
#' @param x \code{lvm}-object
#' @param var the names of the variables that should be removed
#' @param value the names of the variables that should be removed
#' @param expar should the external parameters be also removed
#' @param restaure should the link be kept while removing the linear predictor
#' @param ... argument passed to \code{clean}
#' 
#' @examples
#'
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' kill(m) <- "x7"
#' 
#' rm <- reduce(m, rm.exo = FALSE)
#'
#' kill(rm) <- ~x6
#' kill(rm, value = "x1")
#' kill(rm, value = c("x2","x3"))
#' kill(rm, value = ~x4+x5)
#' kill(rm, value = ~LPy)
#'
#'
#' m <- lvm.reduced()
#' regression(m) <- y1 ~ x1+x2+x3
#' regression(m) <- y2 ~ z
#' covariance(m) <- y1 ~ y2
#' cancel(m) <- y1~x1+x2+x3
#' m
#' cancel(m) <- "y1 ~~ y2"
#' coef(m)
#' 
#' # see test/testthat/test-cancel.R
# }}}

# {{{ reduce.remove.hook
#' @rdname kill.lp
#' @export
lavaReduce.remove.hook  <- function(x, var, expar = TRUE, restaure = FALSE,
                                    ...){
  
  ## lvm object with linear predictor
  if("lvm.reduced" %in% class(x)){
    
    ## remove linear predictors
    if(any(var %in% lp(x))){
      lp.rm <- var[var %in% lp(x)]
      x <- cancel(x, value = lp(x, lp = lp.rm, type = "link"),
                  expar = expar, restaure = restaure, clean = FALSE)    
    }
    
    ## remove variables in the linear predictor
    if(any(var %in% lp(x, type = "x"))){ 
      var.rm <- var[var %in% lp(x, type = "x")]
      index.x <- lapply(lp(x, type = "x", format = "list2"), function(var){which(var$x %in% var.rm)})
      name.lp <- lp(x, type = "name")
      
      for(iterLP in 1:length(index.x)){
        if(length(index.x)==0){next}
        x <- cancel(x, value = lp(x, lp = name.lp[iterLP], type = "link")[index.x[[iterLP]]],
                    expar = expar, restaure = restaure, clean = FALSE)
      }
      
    }
    
  }
  
  return(x)
}

# }}}

# {{{ lavaReduce.cancel.hook
#' @name kill.lp
#' @export
lavaReduce.cancel.hook  <- function(x, value,
                                    expar = TRUE, restaure = FALSE, ...){
  
   if("lvm.reduced" %in% class(x)){
    
    ## normalize value        
    if(is.character(value) && all(value %in% vars(x,lp=TRUE,xlp=TRUE))){
      ## arguments coming from cancel.lvm
      allCoef <- initVar_link(value[1],value[2], format = "txt.formula")
    }else{
      ## argument coming from lavaReduce.remove.hook
      allCoef <- initVar_links(value, format = "txt.formula")
    }
    
    ## coefficients in LP
    allCoef.lp <- intersect(allCoef, lp(x, "link"))
    if(length(allCoef.lp)>0){
      linkLP <- lp(x, type = "link", format = "list")        
      name.lp <- lp(x, type = "name")        
      
      for(iterLP in 1:length(linkLP)){ # iterLP <- 1
        iCoef.lp <- allCoef.lp[allCoef.lp %in% linkLP[[iterLP]]]
        if(length(iCoef.lp)==0){next} 
        
        ## remove external parameters from the LVM 
        if(expar){
          parameter(x, remove = TRUE) <- unlist(iCoef.lp) #must be a character
        }
        
        newlp <- lp(x, type = c("x","con","name","link"), lp = name.lp[iterLP], format = "list2")
        newlp.link <- newlp[[1]]$link
        newlp.x <- newlp[[1]]$x
        
        indexRM <- which(newlp.link %in% iCoef.lp)
        
        ## restaure the links removed from the lp in the lvm               
        if(restaure){                   
          f <- as.formula(paste( lp(x, type = "endo", lp = iterLP), "~",  paste(newlp.x[indexRM], collapse = " + ") ) )
          regression(x) <- f
        }
        
        newlp[[1]]$link <- NULL # link is not part of the lp object. It has been constructed by the lp() extractor so it should be removed for updating using lp<-
        newlp[[1]]$con <- newlp[[1]]$con[-indexRM]
        newlp[[1]]$x <- newlp[[1]]$x[-indexRM]
        
        lp(x, lp = name.lp[iterLP]) <- newlp      
      }
    }
  }
  
  return(x)
}
# }}}
