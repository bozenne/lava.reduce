#' @title Extract variables from a reduced latent variable model
#' @name getReduce
#' @description Extract variables as in a standard lvm but including or not the variables that compose the linear predictor
#' 
#' @param x \code{lvm}-object
#' @param lp should the name of the linear predictors be returned?
#' @param xlp should the name of the variables that the linear predictors aggregates be returned?
#' @param type slot to be return. Can be \code{"link"}, \code{"x"}, \code{"con"}, \code{"name"}, 
#' @param index which linear predictor to consider, 
#' @param format should the results be kept as the list or returned as a single vector, 
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
#' lp(m)
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
#' lp(m)
#' lp(m, type = "x", format = "list")
#' lp(m, index = 1)
#' 
#' exogenous(m)
#' exogenous(m, lp = FALSE)
#' 
#' endogenous(m)
#' endogenous(m, lp = FALSE) # should not change
#' @export
`lp` <- function(x,...) UseMethod("lp")

#' @rdname getReduce 
lp.lvm.reduced <- function(x, type = "name", lp = NULL, format = "vector", ...){
  
  if(length(x$lp)==0){return(NULL)} 
  validNames <- c("link","con","name","x") # names(x$lp[[1]])
  
  ## type
  if(is.null(type)){
    type <- validNames
    size <- FALSE
    format <- "list2"
  }else if(length(type) == 1 && type == "endogeneous"){
    return(names(x$lp))
  }else if(length(type) == 1 && type == "n.link"){
    type <- "link"
    size <- TRUE
    format <- "list"
  }else{
    if(type %in% validNames == FALSE){
      stop("type ",type," is not valid \n",
           "valid types: \"",paste(validNames, collapse = "\" \""),"\" \n")
    }
    size <- FALSE
  }
  
  ## add links
  if("link" %in% type){
    
    for(iterLP in names(x$lp)){
      if(length(x$lp[[iterLP]]$x)>0){
        x$lp[[iterLP]]$link <- paste(iterLP,x$lp[[iterLP]]$x,sep=lava.options()$symbol[1])
      }else{
        x$lp[[iterLP]]$link <- NULL
      }
    }
    
  }
  
  ## format
  validFormat <- c("vector","list","list2")
  if(format %in% validFormat == FALSE){
    stop("format ",format," is not valid \n",
         "format must be on of : \"",paste(validFormat, collapse = " ")," \n")
  }
  if(format != "list2" && length(type)>1){
    stop("format must be \"list2\" when length(type) is not one \n",
         "length(type): ",length(type)," \n")
  }
  
  ## select lp
  if(is.null(lp)){
    lp <- seq_len(length(x$lp))
  }else if(is.numeric(lp)){
    vec <- seq_len(length(x$lp))
    if(any(lp %in% vec == FALSE)){
      stop("lp ",paste(lp, collapse = " ")," is not valid \n",
           "if numeric lp must be in: \"",paste(vec, collapse = " ")," \n")
    }
  }else if(is.character(lp)){
    vec <- unlist(lapply(x$lp, function(x)x[["name"]]))
    if(any(lp %in% vec == FALSE)){
      stop("lp ",paste(lp, collapse = " ")," is not valid \n",
           "if character lp must be in: \"",paste(vec, collapse = "\" \""),"\" \n")
    }
    lp <- match(lp, vec)
  }else{
    stop("lp must be a numeric or character vector \n")
  }
  
  ## extract 
  if(format == "list"){
    res <- lapply(x$lp[lp], function(x)x[[type]])
    
    if(size){
      res <- unlist(lapply(res, function(x) length(x)))
      names(res) <- unlist(lapply(x$lp[lp], function(x)x[["name"]]))
    }
    
  }else{
    res <- lapply(x$lp[lp], function(x)x[type])
  }
  
  ## export
  if(format == "vector"){
    res <- unlist(res)
  }
  
  return(res)
  
}

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

#' @rdname getReduce 
#' @export
endogenous.lvm.reduced <- function(x, top = FALSE, latent = FALSE, ...){
  
  observed <- manifest(x, lp = FALSE, xlp = FALSE)
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

#' @rdname getReduce
#' @export
manifest.lvm.reduced <- function(x, lp = TRUE, xlp = FALSE, ...) {
  
  vars <- vars(x, lp = lp, xlp = xlp)
  if (length(vars)>0)
    setdiff(vars,latent(x))
  else
    NULL
}