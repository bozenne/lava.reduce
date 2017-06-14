#' @title Add regression association to a latent variable model with linear predictor
#' @name regression
#' @description Same as regression for latent variable model with an argument reduce to introduce linear predictors
#' 
#' @param object \code{lvm}-object
#' @param to Character vector of outcome(s) or formula object.
#' @param from Character vector of predictor(s).
#' @param y	Alias for 'to'
#' @param x	Alias for 'from'
#' @param reduce should the from variable be grouped into a linear predictor
#' @param value A formula specifying the linear constraints or if \code{to=NULL} a \code{list} of parameter values.
#' @param ... additional arguments to be passed to \code{regression.lvm}
#' 
#' @examples 
#' 
#' m <- lvm.reduced()
#' m <- regression(m,y='y1',x='x'%++%1:2)
#' m <- regression(m,y='y1',x='z'%++%1:5, reduce = TRUE)
#' m
#' 
#' m <- lvm.reduced()
#' regression(m) <- y1 ~ x1 + x2
#' regression(m, reduce = TRUE) <- y1 ~ x5 + x6 + x7
#' regression(m, reduce = TRUE) <- as.formula(paste0("y~",paste("x",1:5,collapse="+",sep="")))
#'
#' m <- lvm.reduced()
#' regression(m, reduce = "LL") <- y1 ~ x5 + x6 + x7
#' m
#' 

#' @export
#' @rdname regression
`regression<-.lvm.reduced` <- function(object = lvm.reduced(), ..., value){
  
  if(identical(reduce,TRUE)  || is.character(reduce)){
    return(regression(object, value = value, ...))
  }else{
    return(callS3methodParent(object, FUN = "regression<-", class = "lvm.reduced", value = value, ...))
  }
  
}

#' @rdname regression
#' @export
regression.lvm.reduced <- function(object = lvm.reduced(), to, from, y, x, reduce = FALSE, value, ...){

    if(identical(reduce,TRUE)  || is.character(reduce)){

    #### extract to and from
    if (!missing(y)) {
      if (inherits(y, "formula")) 
        y <- all.vars(y)
      to <- y
    }
    if (!missing(x)) {
      if (inherits(x, "formula")) 
        x <- all.vars(x)
      from <- x
    }
    
    #### specific case - see regression.lvm
    if (missing(to)) {
      return(regfix(object))
    }
    if (inherits(to, "formula")) {
      if (!missing(value)) {
        lava::regression(object, to, reduce = reduce, ...) <- value
      }
      else {
        lava::regression(object, reduce = reduce, ...) <- to
      }
      object$parpos <- NULL
      return(object)
    }
    if (is.list(to)) {
      for (t in to) lava::regression(object, reduce = reduce, ...) <- t
      object$parpos <- NULL
      return(object)
    }
    
    #### reduce
    if(any(from %in% lava::latent(object))){
      stop("cannot integrate latent variables in the linear predictor \n",
           "latent variable: \"",paste(from[from %in% lava::latent(object)], collapse = "\" \""),"\" \n")
    }
    
    if(any(from %in% lava::endogenous(object))){
      stop("cannot integrate endogenous variables in the linear predictor \n",
           "endogenous: \"",paste(from[from %in% lava::endogenous(object)], collapse = "\" \""),"\" \n")
    }
    if(!identical(reduce,TRUE) && length(reduce)!=length(to)){
      stop("wrong specification of argument \'reduce\' \n",
           "must be TRUE or a character vector of size ",length(to),"\n")
    }
    
      if("lp" %in% names(object) == FALSE){object$lp <- list()}
     
      for(iterR in to){ ### I don't know where to get the constrains
          name.LP <- if(reduce==TRUE){paste0("LP",iterR)}else{reduce}
          lava::regression(object, to = iterR, from = name.LP, reduce = FALSE, ...) <- 1
          if(iterR %in% names(object$lp)){ 
              # object$lp[[iterR]]$indexCoef <- match(object$lp[[iterR]]$indexCoef, coef(object))
              object$lp[[iterR]]$x <- c(object$lp[[iterR]]$x, from)
              object$lp[[iterR]]$con <- c(object$lp[[iterR]]$con, rep(NA,length(from)))
          }else{
              # object$lp[[iterR]]$indexCoef <- match(namesCoef, coef(object))
              object$lp[[iterR]]$x <- from
              object$lp[[iterR]]$con <- rep(NA,length(from))
              object$lp[[iterR]]$name <- name.LP
          }
          allCoef <- paste0(iterR, lava.options()$symbols[1], object$lp[[iterR]]$x)
          lava::parameter(object) <- c(lava::parameter(object), setdiff(allCoef,coef(object)))
          
          
      }
      return(object)
    
  } else {
    return(callS3methodParent(object, FUN = "regression", class = "lvm.reduced", to = to, from = from, y = y, x = x, value = value, ...))
  }
  
}
