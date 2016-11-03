
#' @title Extract the linear predictors
#' @name lp
#' @description Extract the linear predictors of a lvm object
#' 
#' @param x \code{lvm}-object
#' @param type slot to be return. Can be \code{"link"}, \code{"x"}, \code{"con"}, \code{"name"}
#' @param lp which linear predictor to consider
#' @param format should the results be kept as the list or returned as a single vector
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @examples 
#' ## regresssion
#' m <- lvm.reduced()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' lp(m)
#' lp(m, type = c("x","link"), format = "list2")
#' lp(m, type = NULL)
#' 
#' ## lvm
#' m <- lvm.reduced()
#' m <- regression(m, x=paste0("x",1:10),y="y1", reduce = TRUE)
#' m <- regression(m, x=paste0("x",51:150),y="y2", reduce = TRUE)
#' covariance(m) <- y1~y2
#' 
#' lp(m)
#' lp(m, type = "x", format = "list")
#' lp(m, lp = 1, type = "link")
#' @rdname lp 
#' @export
`lp` <- function(x,...) UseMethod("lp")

#' @rdname lp 
#' @export
lp.lvm.reduced <- function(x, type = "name", lp = NULL, format = "vector", ...){
 
  if(length(x$lp)==0){return(NULL)} 
  validNames <- c("link","con","name","x","endo") # names(x$lp[[1]])
  
  ## type
  if(is.null(type)){
    type <- c("con","name","x")
    size <- FALSE
    format <- "list2"
  }else if(length(type) == 1 && type == "endogeneous"){
    return(names(x$lp))
  }else if(length(type) == 1 && type == "n.link"){
    type <- "link"
    size <- TRUE
    format <- "list"
  }else{
    if(any(type %in% validNames == FALSE)){
      stop("type \"",paste(type[type %in% validNames == FALSE], collapse = "\" \""),"\" is not valid \n",
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
  
  if("endo" %in% type){
    
    for(iterLP in names(x$lp)){
        x$lp[[iterLP]]$endo <- iterLP
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


#' @title Update linear predictors
#' @name updateLP
#' @description Update linear predictors of a lvm object
#' 
#' @param x \code{lvm}-object
#' @param lp the name of the linear predictors to be updated
#' @param value the value that will be allocated
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @examples 
#' ## regresssion
#' m <- lvm.reduced()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' 
#' newLP <- lp(m, type = NULL)[[1]]
#' newLP$link <- newLP$link[1:3]
#' newLP$con <- newLP$con[1:3]
#' newLP$x <- newLP$x[1:3]
#' 
#' lp(m, lp = 1) <- newLP
#' lp(m, type = NULL)
#' 
#' 

#' @rdname updateLP
#' @export
`lp<-` <- function(x, ..., value) UseMethod("lp<-")

#' @rdname updateLP 
#' @export
`lp<-.lvm.reduced` <- function(x, lp = NULL, ..., value){
  
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
    
    if(length(lp) == 1 && length(value) == 3 && all(names(value) == c("con","name","x"))){
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
#' @param link,value the name of the links that should be removed
#' @param lp the name of the linear predictors that should be removed
#' @param expar should the external parameters be also removed
#' @param restaure should the link be kept while removing the linear predictor
#' @param simplify should the class \code{lvm.reduced} be removed is from the \code{lvm}-object if it contains no LP.
#' @param ... argument passed to \code{cancelLP.lvm.reduced}
#' 
#' @examples 
#' ## regresssion
#' m <- lvm.reduced()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' vars(m)
#' 
#' cancelLP(m)
#' cancelLP(m, restaure = TRUE)
#' 
#' lp.links <- lp(m, type = "link")
#'  newm <- cancelLP(m, lp = "LPy", link = lp.links[1:2])
#' lp(newm, type = "link")

#' @rdname cancelLP 
#' @export
`cancelLP` <-
  function(x,...) UseMethod("cancelLP")

#' @rdname cancelLP
#' @export
`cancelLP<-` <- function (x, ..., value) {
  UseMethod("cancelLP<-", x)
}

#' @rdname cancelLP
#' @export
`cancelLP<-.lvm.reduced` <- function (x, ..., value) {
  return(cancelLP(x, link = value, ...))
}

#' @rdname cancelLP
#' @export
cancelLP.lvm.reduced  <- function(x, link = NULL, lp = NULL, expar = TRUE, restaure = FALSE, simplify = TRUE, ...){
  
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
    
    newlp <- lp(x, type = c("x","con","name","link"), lp = names.lp[iterLP], format = "list2")
    newlp.link <- newlp[[1]]$link
    newlp.x <- newlp[[1]]$x
    
    indexRM <- which(newlp.link %in% link[[iterLP]])
    
    if(length(indexRM) != length(link[[iterLP]])){
      stop("unknown links: \"",paste(link[[iterLP]][link[[iterLP]] %in% newlp.link == FALSE], collapse = "\" \""),"\" \n",
           "possible links for ",names.lp[iterLP],": \"",paste(newlp.link[newlp.link %in% link[[iterLP]] == FALSE], collapse = "\" \""),"\" \n")
    }
    
    if(restaure){
      f <- as.formula(paste( lp(x, type = "endo", lp = iterLP), "~",  paste(newlp.x[indexRM], collapse = " + ") ) )
      regression(x) <- f
    }


    newlp[[1]]$link <- NULL
    newlp[[1]]$con <- newlp[[1]]$con[-indexRM]
    newlp[[1]]$x <- newlp[[1]]$x[-indexRM]
    
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

