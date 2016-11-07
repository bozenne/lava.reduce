#' @title Remove variables from linear predictors
#' @description Remove one or several variable from the linear predictors
#' @name killLP
#' 
#' @param x \code{lvm}-object
#' @param value the name of the links that should be removed
#' @param lp should the variables be removed from the linear predictor
#' @param expar should the external parameters be also removed
#' @param restaure should the link be kept while removing the linear predictor
#' @param simplify should the class \code{lvm.reduced} be removed is from the \code{lvm}-object if it contains no LP.
#' @param ... argument passed to \code{kill.lvm.reduced}
#' 
#' @examples 
#' # see test/testthat/test-cancel.R

#' @rdname killLP
#' @export
kill.lvm.reduced  <- function(x, value = NULL, lp = TRUE, expar = TRUE, restaure = FALSE, simplify = TRUE, ...){
  
  if("formula" %in% class(value)){
    value <- list(value)
  }
  allVar <- unlist(lapply(value, function(f){unlist(initVar_link(f))}))
  
  if(is.null(lp(x))){
    return(kill.lvm(x, value = value, ...))
  }
  
  if(any(allVar %in% lp(x))){
    lp.rm <- allVar[allVar %in% lp(x)]
    x <- cancel(x, value = lp(x, lp = lp.rm, type = "link"), expar = expar, restaure = restaure, simplify = simplify, ...)
    
  }
  
  if("lvm.reduced" %in% class(x) && lp && any(allVar %in% lp(x, type = "x"))){ ## standard
    var.rm <- allVar[allVar %in% lp(x, type = "x")]
    index.x <- lapply(lp(x, type = "x", format = "list2"), function(var){which(var$x %in% var.rm)})
    name.lp <- lp(x, type = "name")
    
    for(iterLP in 1:length(index.x)){
      if(length(index.x)==0){next}
      x <- cancel(x, value = lp(x, lp = name.lp[iterLP], type = "link")[index.x[[iterLP]]], expar = expar, restaure = restaure, simplify = simplify, ...)
    }
    
  }
  
  if(any(allVar %in% vars(x, lp = FALSE, xlp = FALSE))){ ## standard
    var.rm <- allVar[allVar %in% vars(x, lp = FALSE, xlp = FALSE)]
    
    x <- kill.lvm(x, value = var.rm, ...)
    
  }
  
  return(x)
}


#' @name killLP
#' @export
cancel.lvm.reduced  <- function(x, value = NULL, expar = TRUE, restaure = FALSE, simplify = TRUE, ...){
  
  if("formula" %in% class(value)){
    value <- list(value)
  }
  
  ## order by linear predictor
  allCoef <- lapply(value, function(link){initVar_link(link, format = "txt.formula")})
  
  if(any(unlist(allCoef) %in% coef(x) == FALSE)){
    stop("unknown link: \"",paste(allCoef[allCoef %in% coef(x) == FALSE], collapse = "\" \""),"\" \n",
         "possible link: \"",paste(coef(x)[coef(x) %in%allCoef == FALSE], collapse = "\" \""),"\" \n")
  }
  
  allCoef.lp <- intersect(unlist(allCoef), lp(x, "link"))
  allCoef.nlp <- setdiff(unlist(allCoef), lp(x, "link"))
  
  ## normal coefficients
  if(length(allCoef.nlp)>0){
    
    for(iterF in 1:length(allCoef.nlp)){
      
      x <- cancel.lvm(x, value = unlist(initVar_link(allCoef.nlp[[iterF]])))
    }
  }
  
  ## coefficients in LP
  if(length(allCoef.lp)>0){
    
    vecLP.response <- sapply(allCoef.lp, function(link){select.response(character2formula(link))})
    value <- tapply(allCoef.lp, vecLP.response, list)
    
    names.lp <- lp(x, type = "name")[match(unique(vecLP.response), lp(x, type = "endo"))]
    n.lp <- length(names.lp)
    
    for(iterLP in 1:n.lp){
      
      ## remove external parameters from the LVM 
      if(expar){
        parameter(x, remove = TRUE) <- unlist(value[[iterLP]]) #must be a character
        # x$expar <- x$expar[setdiff(names(x$expar),coefLP)]
        # x$exfix <- x$exfix[setdiff(names(x$exfix),coefLP)]
        # if(length(x$exfix)==0){x$exfix <- NULL}
        # x$attributes$parameter <- x$attributes$parameter[setdiff(names(x$attributes$parameter),coefLP)]
        # x$index$npar.ex <- x$index$npar.ex[setdiff(names(x$index$npar.ex),coefLP)]
        # x$index$parval <- x$index$parval[setdiff(names(x$index$parval),coefLP)]
      }
      
      newlp <- lp(x, type = c("x","con","name","link"), lp = names.lp[iterLP], format = "list2")
      newlp.link <- newlp[[1]]$link
      newlp.x <- newlp[[1]]$x
      
      indexRM <- which(newlp.link %in% value[[iterLP]])
      
      if(restaure){
        f <- as.formula(paste( lp(x, type = "endo", lp = iterLP), "~",  paste(newlp.x[indexRM], collapse = " + ") ) )
        regression(x) <- f
      }
      
      
      newlp[[1]]$link <- NULL
      newlp[[1]]$con <- newlp[[1]]$con[-indexRM]
      newlp[[1]]$x <- newlp[[1]]$x[-indexRM]
      
      lp(x, lp = names.lp[iterLP]) <- newlp
      
    }
  }
  
  ## remove empty LP
  n.link <- lp(x, type = "n.link")
  if(any(n.link==0)){
    name.lp <- names(n.link)[which(n.link==0)]
    x <- kill.lvm(x, name.lp)
    x$lp[lp(x, lp = name.lp, type = "endo")] <- NULL
  }
  
  ## update class
  if(simplify && length(lp(x, type = "link", format = "vector")) == 0 ){
    class(x) <- setdiff(class(x), "lvm.reduced")
  }
  
  return(x)
}

