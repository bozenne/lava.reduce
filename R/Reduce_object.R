#### intialisation ####

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

#### convertion ####

#' @title Conversion to a latent variable model with linear predictors
#' @description Add the possibility to use linear predictor in a latent variable model
#' 
#' @param x \code{lvm}-object
#' 
lvm2reduce <- function(x){
  
  if("lvm.reduced" %in% class(x) == FALSE){
    
    x$lp <- initLP()
    class(x) <- append("lvm.reduced",class(x))
    
  }else{
    warning("x is already a reduced latent variable model \n")
  }
  
  return(x)
}

#### specification ####

#' @title Reduce latent variable model
#' @name reduce
#' @description Aggregate exogeneous variables into a linear predictor
#' 
#' @param x \code{lvm}-object
#' @param link the links that should be aggregated into linear predictors
#' @param endo the endogeneous variables for which the related exogeneous variables should be aggregated
#' @param rm.exo should the exogeneous variables be remove from the object
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @examples 
#' ## regression
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' 
#' m
#' reduce(m)
#' coef(reduce(m))
#' reduce(m, rm.exo = FALSE)
#' 
#' ## lvm
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y1")
#' m <- regression(m, x=paste0("x",11:30),y="y2")
#' covariance(m) <- y1~y2
#' 
#' m
#' reduce(m)
#' coef(reduce(m))
#' reduce(m, rm.exo = FALSE)
#'
#' @export
`reduce` <-
  function(x,...) UseMethod("reduce")

#' @rdname reduce
#' @export
reduce.lvm <- function(x, link = NULL, endo = NULL, rm.exo = TRUE, ...){
  
  if("lvm.reduced" %in% class(x) == FALSE){
    x <- lvm2reduce(x)
  }
  
  #### define the links
  if(!is.null(link)){ # reduce specific links
    
    if("formula" %in% class(link)){ link <- list(link) }
    ls.link <- lapply(link, initVar_link, repVar1 = TRUE)
    vec.endo <- unlist(lapply(ls.link, "[[", 1))
    vec.exo <- unlist(lapply(ls.link, "[[", 2))
    
    exo <- tapply(vec.exo, vec.endo, list)
    endo <- names(exo)
   
  }else{ # reduce the linear predictor of specific endogeneous variables
    if(is.null(endo)){endo <- lava::endogenous(x)}
    
    # if nothing to reduce return object
    remaining.exo <- intersect(vars(x, lp = FALSE), exogenous(x, lp = FALSE, xlp = TRUE))
    if(length(remaining.exo)==0){return(x)}
    
    M <- x$M[remaining.exo,endo, drop = FALSE]
    
    col.reg <- names(which(apply(M, 2, function(i){any(i==1)}))) # find endo having exo
    exo <- lapply(col.reg, function(var){rownames(M)[which(M[,var,drop = FALSE]==1)]})
    names(exo) <- col.reg
    
    if(is.null(names(exo))){
      cat("no regression model has been found in the object \n")
      return(x)
    }else{
      endo <- names(exo)#vars(x)[index.reduce]
    }
    
  }
  
  #### reduce
  n.endo <- length(endo)
  
  for(iterR in 1:n.endo){
    name.endo <- endo[iterR]
    name.exo <- exo[[iterR]]
    
    if(length(name.exo)>0){
      ## can be problematic as we don't know about "additive" or other possibly relevant arguments
      f <- stats::as.formula(paste(name.endo,"~",paste(name.exo, collapse = "+")))
      cancel(x) <- f
      
      x <- lava::regression(x, to = name.endo, from = name.exo, reduce = TRUE)
    }
  }
  
  if(rm.exo){
    indexClean <- which(rowSums(x$index$A[x$exogenous,,drop = FALSE]!=0)==0)
    kill(x) <- x$exogenous[indexClean]
  }
  
  return(x)
}
