# scoreMVN_mets <- function(Y, mu, dmu, S, dS){
#   
#   .Call("loglikMVN", yl = Y, yu = NULL, status = NULL, 
#         mu = mu, dmu = dmu, s = S, 
#         ds = dS, z = NULL, su = NULL, dsu = NULL, 
#         threshold = NULL, dthreshold = NULL, score = TRUE, PACKAGE = "mets")
#   
# }

#' @title Moments for a reduced lvm model
#' @name momentLVMr
#' @description Compute the likelihood, score and hessian of a reduced lvm model
#' 
#' @param x,object \code{lvm.reduced}-object
#' @param p current parameter estimate
#' @param data dataset
#' @param n number of observations
#' @param indiv should the individual contribution be returned? Else average over the observations
#' @param type method used to compute the hessian
#' @param implementation default is R. If set to cpp all the computation of the gradient is made in a C++ routine.
#' @param ... additional arguments
#' 
#' @details this function assumes that the external parameters in p are at the end of the vector
#' 

#' @rdname momentLVMr
#' @export
gaussianLP_method.lvm <- "nlminb2"

#' @rdname momentLVMr
#' @export
gaussianLP_objective.lvm <- function(x, p, data, ...){ 
  l <- gaussianLP_logLik.lvm(x, p=p, data=data, ...)
  return(-l)
}


#' @rdname momentLVMr
#' @export
gaussianLP_logLik.lvm <- function(object, p, data, ...)  {
  
  if(is.matrix(p)){p <- as.double(p)}
  if(is.null(names(p))){names(p) <- coef(object)}
  
  dataLP <- calcLP(data, p = p, 
                   lp.x = lp(object, type = "x", format = "list"), 
                   lp.link = lp(object, type = "link", format = "list"), 
                   lp.name = lp(object, type = "name"))
  
  ## from normal_objective.lvm
  y <- lava::index(object)$endogenous
  ord <- lava::ordinal(object)
  status <- rep(0,length(y))
  bin <- tryCatch(match(do.call("binary",list(x=object)),y),error=function(x) NULL)
  status[match(ord,y)] <- 2
  
  mu <- stats::predict(object,data=dataLP,p=p) ## can be optimized
  S <- attributes(mu)$cond.var
  class(mu) <- "matrix"
  thres <- matrix(0,nrow=length(y),max(1,attributes(ord)$K-1)); rownames(thres) <- y
  for (i in seq_len(length(attributes(ord)$fix))) {
    nn <- names(attributes(ord)$idx)[i]
    ii <- attributes(ord)$idx[[nn]]
    val <- (attributes(mu)$e[ii])
    thres[nn,seq_len(length(val))] <-
      cumsum(c(val[1],exp(val[-1])))
  }
  
  yl <- yu <- as.matrix(data[,y,drop=FALSE])
  
  if (!inherits(yl[1,1],c("numeric","integer","logical")) ||
      !inherits(yu[1,1],c("numeric","integer","logical")))
    stop("Unexpected data (normal_objective)")
  
  l <- sum(mets::loglikMVN(yl = yl, yu = yu, status = status, mu = mu, S = S, thres = thres))
  
  return(l)
}

#' @rdname momentLVMr
#' @export
gaussianLP_gradient.lvm <- function(x, p, data, ...){
  
  val <- -gaussianLP_score.lvm(x, p = p, data = data, ...)
  if (!is.null(nrow(val))){
    val <- colSums(val)
  }
  
  return(val)
}

#' @rdname momentLVMr
#' @export
gaussianLP_score.lvm <- function(x, p, data, indiv = FALSE, implementation = "R", ...)  {
  
  
  if(is.matrix(p)){p <- as.double(p)}
  if(is.null(names(p))){names(p) <- coef(x)}
  
  #### extract information
  lp.endo <- lp(x, type = "endogeneous")
  lp.link <- lp(x, type = "link", format = "list")
  lp.x <- lp(x, type = "x", format = "list")
  
  names.data <- colnames(data)
  names.p <- names(p)
 
  M <- moments(x,p) # lava:::moments.lvm
  D <- deriv.lvm(x,p=p) # [WARNING lava:::]
  
  
  if(implementation == "cpp"){
    
    score <- scoreLVM(data = as.matrix(data),
                      p = p,
                      mu = M$xi%x%rep(1,NROW(data)),
                      S = M$C,
                      dmu = D$dxi,
                      dS = D$dS,
                      scoreFun = mets::scoreMVN,#scoreMVN_mets,
                      indexCoef = lapply(lp.link, function(link){match(link, names.p) - 1}),
                      indexEndo = lapply(lp.x, function(link){match(link, names.data) - 1}),
                      indexIntercept = match(lp.endo, names.p) - 1,
                      indexLP = match(lp(x, type = "name"), names.data) - 1,
                      indexManifest = match(lava::manifest(x), names.data) - 1,
                      indiv = indiv
    )
    colnames(score) <- coef(x)
    
  }else{
    
    #### compute
    dataLP <- calcLP(data, p = p,
                     lp.x = lp.x, lp.link = lp.link, lp.name = lp(x, type = "name"))
    
    score <- mets::scoreMVN(y = dataLP[,lava::manifest(x),drop = FALSE],
                            mu = M$xi%x%rep(1,NROW(dataLP)),
                            S= M$C,
                            dmu = D$dxi,
                            dS = D$dS)
    colnames(score) <- coef(x)#c(name.intercept,name.regression,name.covariance)
   
    
    ## apply chain rule
    n.lp <- length(x$lp)
    
    for(iterLP in 1:n.lp){ 
      name.intercept <- lp.endo[iterLP]
      
      ## extract data
      X <-  data[,lp.x[[iterLP]],drop = FALSE]
      dlp <- X # - what about dB/db in presence of constrains
      score[,lp.link[[iterLP]]] <- apply(dlp, 2, function(j){j*score[,name.intercept]})
    }
    
    if(indiv == FALSE){
      score <- rbind(apply(score,2,sum))
    }
  }
  
  return(score) 
}

#' @rdname momentLVMr
#' @export
gaussianLP_hessian.lvm <- function(x, p, n, type,...) {
  dots <- list(...);
  
  if(type == "E"){
    S <- -gaussianLP_score.lvm(x,p=p,data=dots$data,indiv = TRUE)
    I <- t(S)%*%S
    attributes(I)$grad <- colSums(S)
    
    return(I)
  }else if(type=="num"){
    myg <- function(p1){ -gaussianLP_score.lvm(x,p=p1,n=n,data=dots$data,indiv=FALSE) }
    I <- numDeriv::jacobian(myg,p) 
    I <- (I+t(I))/2
    return( I )
  }else if(type == "information"){ ## true part
    
    lp.name <- lp(x, type = "name")
    lp.link <- lp(x, type = "link", format = "list")
    lp.x <- lp(x, type = "x", format = "list")
    
    dataLP <- calcLP(dots$data,  p = p,
                     lp.x = lp.x, lp.link = lp.link, lp.name = lp.name)
    
    ## direct
    I <- information(x=x, p=p, n=n, data = dataLP)
    return(I)
    
  }
  
}

#' @rdname momentLVMr
#' @export
gaussian1LP_method.lvm <- gaussianLP_method.lvm
#' @rdname momentLVMr
#' @export
gaussian1LP_logLik.lvm <- gaussianLP_logLik.lvm
#' @rdname momentLVMr
#' @export
gaussian1LP_objective.lvm <- gaussianLP_objective.lvm
#' @rdname momentLVMr
#' @export
gaussian1LP_score.lvm <- gaussianLP_score.lvm
#' @rdname momentLVMr
#' @export
gaussian1LP_gradient.lvm <- gaussianLP_gradient.lvm
#' @rdname momentLVMr
#' @export
gaussian1LP_hessian.lvm <- function(x, type, ...){
  gaussianLP_hessian.lvm(x, type = "num", ...)
}

#' @rdname momentLVMr
#' @export
gaussian2LP_method.lvm <- gaussianLP_method.lvm
#' @rdname momentLVMr
#' @export
gaussian2LP_logLik.lvm <- gaussianLP_logLik.lvm
#' @rdname momentLVMr
#' @export
gaussian2LP_objective.lvm <- gaussianLP_objective.lvm
#' @rdname momentLVMr
#' @export
gaussian2LP_score.lvm <- gaussianLP_score.lvm
#' @rdname momentLVMr
#' @export
gaussian2LP_gradient.lvm <- gaussianLP_gradient.lvm
#' @rdname momentLVMr
#' @export
gaussian2LP_hessian.lvm <- function(x, type, ...){
  gaussianLP_hessian.lvm(x, type = "E", ...)
}



#' @title Compute the linear predictor
#' @description Compute the value of the linear predictors of a LVM and store it into the dataset
#' 
#' @inheritParams momentLVMr
#' @param lp.x the name of the exogeneous variables involved in the LP
#' @param lp.link list containing for each LP the name of the links
#' @param lp.name list containing for each LP its name
calcLP <- function(data, p,
                   lp.x, lp.link, lp.name){
  
  n.lp <- length(lp.name)
  
  for(iterLP in 1:n.lp){
    
    ## extract coefficients according to constrains
    b <- p[lp.link[[iterLP]]] # x$lp$y1$con
    
    ## compute linear predictor for the reduce model
    X <- data[,lp.x[[iterLP]],drop = FALSE]
    data[,lp.name[iterLP]] <- as.matrix(X) %*% b
    
  }
  return(data)
}
