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
#' @param endo the endogeneous variables for which the related exogeneous variables should be aggregated
#' @param rm.exo should the exogeneous variables be remove from the object
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

#' @rdname getReduce
reduce.lvm <- function(object, link = NULL, endo = NULL, rm.exo = TRUE){
  
  if("lvm.reduced" %in% class(object) == FALSE){
    object <- lvm2reduce(object)
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
    if(is.null(endo)){endo <- endogenous(object)}
    
    # if nothing to reduce return object
    remaining.exo <- intersect(vars(object, lp = FALSE), exogenous(object, lp = FALSE, xlp = TRUE))
    if(length(remaining.exo)==0){return(object)}
    
    M <- object$M[remaining.exo,endo, drop = FALSE]
    
    col.reg <- names(which(apply(M, 2, function(x){any(x==1)}))) # find endo having exo
    exo <- lapply(col.reg, function(var){names(which(M[,var]==1))})
    names(exo) <- col.reg
    
    if(is.null(names(exo))){
      cat("no regression model has been found in the object \n")
      return(object)
    }else{
      endo <- names(exo)#vars(object)[index.reduce]
    }
    
    
  }
  
  #### reduce
  n.endo <- length(endo)
  
  for(iterR in 1:n.endo){
    name.endo <- endo[iterR]
    name.exo <- exo[[iterR]]
    
    if(length(name.exo)>0){
      ## can be problematic as we don't know about "additive" or other possibly relevant arguments
      f <- as.formula(paste(name.endo,"~",paste(name.exo, collapse = "+")))
      cancel(object) <- f
      
      object <- regression(object, to = name.endo, from = name.exo, reduce = TRUE)
    }
  }
  
  if(rm.exo){
    indexClean <- which(rowSums(object$index$A[object$exogenous,]!=0)==0)
    kill(object) <- object$exogenous[indexClean]
  }
  
  return(object)
}

#### update/cancel ####

#' @title Update of the linear predictor
#' @name lp
#' @description Update one or several linear predictors
#' 
#' @param x \code{lvm}-object
#' @param lp the name of the linear predictors to be updated
#' @param value the value that will be allocated
#' 
#' @examples #' 
#' ## regresssion
#' m <- lvm.reduced()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' vars(m)
#' 
#' newLP <- lp(m, type = NULL)[[1]]
#' newLP$link <- newLP$link[1:3]
#' newLP$con <- newLP$con[1:3]
#' newLP$x <- newLP$x[1:3]
#' 
#' lp(m, lp = 1) <- newLP
#' lp(m, type = NULL)
#' 
#' @export
`lp<-` <- function(x,...) UseMethod("lp<-")

#' @rdname lp 
#' @export
`lp<-.lvm.reduced` <- function(x, lp = NULL, value){
  
  if(is.null(lp)){
    
    x$lp <- value
    
  }else{
    
    ## valid LP
    if(is.numeric(lp)){
      vec <- seq_len(length(x$lp))
      if(any(lp %in% vec == FALSE)){
        stop("lp ",paste(lp, collapse = " ")," is not valid \n",
             "if numeric lp must be in: \"",paste(vec, collapse = " ")," \n")
      }
    }else if(is.character(lp)){
      vec <- unlist(lapply(x$lp, function(x)x[["name"]]))
      if(any(lp %in% vec == FALSE)){
        stop("lp ",paste(lp, collapse = " ")," is not valid \n",
             "if character lp must be in: \"",paste(vec, collapse = " ")," \n")
      }
      lp <- match(lp, vec)
    }
    
    if(length(lp) == 1 && length(value) == 4 && all(names(value) == c("link","con","name","x"))){
      value <- list(value)
    }
    
    x$lp[lp] <- value
    
  }
  
  return(x)
  
}

#' @title Remove variable from the linear predictor
#' @description Remove one or several variable from the linear predictor
#' @name cancelLP
#' 
#' @param x \code{lvm}-object
#' @param lp should the name of the variables corresponding to the linear predictors be returned?
#' @param expar should the external parameters be also removed
#' 
#' @examples 

#' ## regresssion
#' m <- lvm.reduced()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' vars(m)
#' 
#' cancelLP(m)
#' lp.links <- lp(m, type = "link")
#'  newm <- cancelLP(m, lp = "LPy", link = lp.links[1:2])
#' lp(newm, type = "link")
#' 
#' @export
`cancelLP` <-
  function(x,...) UseMethod("cancelLP")

#' @export
#' @rdname cancelLP
`cancelLP<-` <- function (x, ..., value) {
  UseMethod("cancelLP<-", x)
}

#' @export
#' @rdname cancelLP
`cancelLP<-.lvm.reduced` <- function (x, ..., value) {
  return(cancelLP(x, link = value, ...))
}

#' @export
#' @rdname cancelLP
cancelLP.lvm.reduced  <- function(x, link = NULL, lp = NULL, expar = TRUE, simplify = TRUE){
  
  if(is.null(lp)){
    lp <- lp(x) 
  }else if(!is.list(lp)){
    lp <- as.list(lp)
  }
  n.lp <- length(lp)
  names.lp <- unlist(lp)
  
  if(is.null(link)){
    link <- lp(x, type = "link", lp = unlist(lp), format = "list")
  }else if(!is.list(link)){
    link <- lapply(1:n.lp, function(x) link)
    names(link) <- names.lp
  }
  
  ## remove external parameters from the LVM 
  if(expar){
    parameter(x, remove = TRUE) <- unlist(link)
    # x$expar <- x$expar[setdiff(names(x$expar),coefLP)]
    # x$exfix <- x$exfix[setdiff(names(x$exfix),coefLP)]
    # if(length(x$exfix)==0){x$exfix <- NULL}
    # x$attributes$parameter <- x$attributes$parameter[setdiff(names(x$attributes$parameter),coefLP)]
    # x$index$npar.ex <- x$index$npar.ex[setdiff(names(x$index$npar.ex),coefLP)]
    # x$index$parval <- x$index$parval[setdiff(names(x$index$parval),coefLP)]
  }
  
  ## removing variable in the LP
  for(iterLP in 1:n.lp){
    newlp <- lp(x, type = NULL, lp = names.lp[iterLP])[[1]]
    
    indexRM <- which(newlp$link %in% link[[iterLP]])
    
    if(length(indexRM) != length(link[[iterLP]])){
      stop("unknown links: \"",paste(link[[iterLP]][link[[iterLP]] %in% newlp$link == FALSE], collapse = "\" \""),"\" \n",
           "possible links for ",names.lp[iterLP],": \"",paste(newlp$link[newlp$link %in% link[[iterLP]] == FALSE], collapse = "\" \""),"\" \n")
    }
    
    newlp$link <- newlp$link[-indexRM]
    newlp$con <- newlp$con[-indexRM]
    newlp$x <- newlp$x[-indexRM]
    
    lp(x, lp = names.lp[iterLP]) <- newlp
  }
 
  ## remove empty LP
  n.link <- lp(x, type = "n.link")
  if(any(n.link==0)){
    rmvar(x) <- names(n.link)[which(n.link==0)]
  }
  
  ## update class
  if(simplify && length(lp(x, type = "link", format = "vector")) == 0 ){
    class(x) <- setdiff(class(x), "lvm.reduced")
  }
  
  return(x)
}

