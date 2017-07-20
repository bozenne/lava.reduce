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

#' @title Remove all linear predictors from a lvm.reduce object
#' @description Remove all linear predictors from a lvm.reduce object. Restaure the links with the covariates.
#' @name reduce2lvm
#' 
#' @param x \code{lvm}-object
#' @param ... not used.
#' 
#' @examples 
#' m <- lvm(Y~X1+X2+X3)
#' mR <- reduce(m)
#' 
#' reduce2lvm(mR)
#' 
#' m <- lvm(c(Y1~X1+X2+X3,Y2~Z1+X3))
#' mR <- reduce(m)
#' regression(mR) <- Y3~X1+X8
#' reduce2lvm(mR)
#'  
#' @export
`reduce2lvm` <-
  function(x,...) UseMethod("reduce2lvm")
# }}}

#' @rdname reduce2lvm
#' @export
reduce2lvm.lvm.reduced <- function(x, ...){
  
  if("lvm.reduced" %in% class(x)){
     
    allLinks <- paste0(lp(x,type = "endo"),lava.options()$symbols[1],lp(x,type = "name"))
    
    for(iLP in 1:length(allLinks)){ # iLP <- 1
      allLinks[iLP] <- gsub(lava.options()$symbols[1],"~",allLinks[iLP])
      cancel(x, restaure = TRUE) <- as.formula(allLinks[iLP])
    }
    
    x <- clean(x, rm.endo = FALSE)
    
  }else{
    warning("x is already a reduced latent variable model \n")
  }
  
  return(x)
}

#' @rdname reduce2lvm
#' @export
reduce2lvm.lvmfit.reduced <- function(x, ...){
   
  data <- model.frame(x)
  xfull <- reduce2lvm(Model(x))
    
  ## from lava::estimate
  resProc <- procdata.lvm(xfull,data=data)
  xfix <- setdiff(colnames(data)[(colnames(data) %in% parlabels(xfull, exo = TRUE))], latent(xfull))
  fix <- ifelse(length(xfix) > 0, FALSE, TRUE)
  xfull <- fixsome(xfull, measurement.fix = fix, S = resProc$S, mu = resProc$mu, 
                   n = resProc$n, debug = FALSE)
  
  if (length(xfix) > 0) {
    index(xfull) <- reindex(xfull, sparse = FALSE, zeroones = TRUE, 
                        deriv = TRUE)
  }
  else {
    xfull <- updatelvm(xfull, sparse = FALSE, zeroones = TRUE, 
                   deriv = TRUE, mean = TRUE)
  }
  
  ## export
  return(xfull)
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
  
  myhooks <- lavaReduce::gethook_lavaReduce("reduce.hooks")
  save <- NULL
  for (f in myhooks) {
    resHook <- do.call(f, list(x=x, link = link, ...))
    if("x" %in% names(resHook)){x <- resHook$x}
    if("link" %in% names(resHook)){link <- resHook$link}
    if("save" %in% names(resHook)){save <- c(save,resHook$save)}
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
    
    ## define all links
    link <- paste(colnames(M)[oneSearch[,2]],
                  rownames(M)[oneSearch[,1]],
                  sep = lava.options()$symbols[1])
    
    ## define exogeneous variables relative to each endogenous variable
    all.y <- colnames(M)
    n.y <- length(all.y)
    ls.x <- lapply(1:n.y, function(iY){
      if(any(oneSearch[,2]==iY)){
        rownames(M)[oneSearch[oneSearch[,2]==iY,1]]
      }
    })
    names(ls.x) <- all.y
    ls.x <- ls.x[lengths(ls.x)>0]
  }else{
    ls.xy <- initVar_links(link, format = "list")
    ls.x <- tapply(ls.xy$var2, ls.xy$var1, function(x){x})
  }
  
  all.y <- names(ls.x)
  n.y <- length(ls.x)
  
 #### define the links
  # can be problematic as we don't know about "additive" or other possibly relevant arguments
  ls.link <- combine.formula(link)
  
  #### reduce
  n.endo <- length(ls.link)
  if("lvm.reduced" %in% class(x) == FALSE){
    x <- lvm2reduce(x)
  }
  
  x.class <- class(x)
  class(x) <- append("lvm.reduced",x.class)    
  
  for(iterR in 1:n.endo){# iterR <- 1        
    cancel(x) <- ls.link[[iterR]]
  }
  if(clean){
    x <- clean(x, rm.endo = FALSE)
  }
  for(iterR in 1:n.endo){# iterR <- 1        
    x <- regression.lvm.reduced(x, reduce = TRUE, y = all.y[iterR], x = ls.x[[all.y[iterR]]])
  }
  class(x) <- x.class
  
  # restaure
  if(!is.null(save)){
    for(iElement in names(save)){
    x[[iElement]] <- save[[iElement]]
    }
  }
  
  return(x)
}
# }}}

