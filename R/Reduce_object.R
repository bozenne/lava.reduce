# {{{ intialisation

#' @title Initialize a linear predictor
#' @description Initialise a linear predictor
#' 
initLP <- function(){
  
  lp <- list()
  
  class(lp) <- "lp"
  return(lp)
  
}

#' @title Initialize new latent variable model with linear predictors
#' @description Function that constructs a new latent variable model object with linear predictors
#' 
#' @param ... arguments to be passed to the lvm function
#' 
#' @examples 
#' m <- lvm.reduced(Y ~ X1+X2)
#' class(m)
#' lp(m)
#' 
#' @export
lvm.reduced <- function(...){
  x <- lvm(...)
  x <- lvm2reduce(x)
  return(x)
}

# }}}

# {{{ conversion
#' @title Conversion to a latent variable model with linear predictors
#' @description Add the possibility to use linear predictor in a latent variable model
#' 
#' @param x \code{lvm}-object
#'
#' @export
lvm2reduce <- function(x){
  
  if("lvm.reduced" %in% class(x) == FALSE){
    
    x$lp <- initLP()
    class(x) <- append("lvm.reduced",class(x))
    
  }else{
    warning("x is already a reduced latent variable model \n")
  }
  
  return(x)
}

# }}}
# {{{ reduce

# {{{ doc
#' @title Reduce latent variable model
#' @name reduce
#' @description Aggregate exogeneous variables into a linear predictor
#' 
#' @param x \code{lvm}-object
#' @param link the links that should be aggregated into linear predictors
#' @param endo the endogeneous variables for which the related exogeneous variables should be aggregated
#' @param clean should the lvm object be simplified using the \code{clean} function
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @examples 
#' ## regression
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' m
#'
#' mtest <- m
#' cancel(mtest) <- y ~ x1 + x2 + x3 + x4 + x5
#' mtest
#' coef(mtest)
#' rm <- reduce(m)
#' rm
#' lp(rm, type = "link")
#' coef(rm)
#' 
#' coef(reduce(m))
#' reduce(m, rm.exo = TRUE)
#' 
#' reduce(m, link = paste0("y~x",1:5))
#' 
#' ## lvm
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y1")
#' m <- regression(m, x=paste0("x",11:30),y="y2")
#' covariance(m) <- y1~y2
#'
#' system.time(
#' rm <- reduce(m)
#' )
#' coef(rm)
#' reduce(m, rm.exo = TRUE)
#'
#' @export
`reduce` <-
  function(x,...) UseMethod("reduce")
# }}}

# {{{ reduce.lvm 
#' @rdname reduce
#' @export
reduce.lvm <- function(x, link = NULL, endo = NULL, clean = TRUE, ...){

   myhooks <- gethook_lava.reduce("reduce.hooks")
   for (f in myhooks) {
       res <- do.call(f, list(x=x, link = link, ...))
       if("x" %in% names(res)){x <- res$x}
       if("link" %in% names(res)){link <- res$link}
   }
    
  if(is.null(link)){ # reduce the linear predictor of specific endogeneous variables
  
    if(is.null(endo)){endo <- lava::endogenous(x)}

    # if nothing to reduce return object
    remaining.exo <- intersect(vars(x, lp = FALSE), 
                               exogenous(x, lp = FALSE, xlp = TRUE))
    if(length(remaining.exo)==0){
      cat("All exogenous variables has already been reduces \n")
      return(x)
    }
    
    M <- x$M[remaining.exo,endo, drop = FALSE]
    oneSearch <- which(M==1, arr.ind = TRUE)
    
    if(NROW(oneSearch) == 0){
      cat("No regression model has been found in the object \n")
      return(x)
    }
    
    link <- paste(colnames(M)[oneSearch[,2]],
                  rownames(M)[oneSearch[,1]],
                  sep = lava.options()$symbols[1])
    
  }
  
    #### define the links
    # can be problematic as we don't know about "additive" or other possibly relevant arguments
    ls.link <- combine.formula(link)
  
    #### reduce
    n.endo <- length(ls.link)
    if("lvm.reduced" %in% class(x) == FALSE){
        x <- lvm2reduce(x)
    }

    x.class <- class(x)
    class(x) <- c("lvm.reduced","lvm")    
    for(iterR in 1:n.endo){# iterR <- 1        
        cancel(x) <- ls.link[[iterR]]
        regression(x, reduce = TRUE) <- ls.link[[iterR]]
    }
    class(x) <- x.class

    if(clean){
        x <- clean(x, ...)
    }
  
  return(x)
}
# }}}

# }}}
