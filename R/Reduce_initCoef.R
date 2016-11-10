#' @title Initialise \code{lvm.reduced}
#' @name initCoef
#' @description Estimate initial values for the parameters contained in a \code{lvm.reduced} object
#' 
#' @inheritParams estimate.lvm.reduced
#' @param optim control/optimization parameters
#' @param ... additional arguments to be passed to the low level functions
#' 

#' @rdname initCoef
#' @export
`initCoef` <-  function(x,...) UseMethod("initCoef")

#' @rdname initCoef
#' @export
initCoef <- function(x, data, optim, ...){
  
  ## compute moment
  dd <- procdata.lvm(x, data = data) # callS3methodParent(x, FUN = "procdata", class = "lvm.reduced", data = data)
  S <- dd$S
  mu <- dd$mu
  n <- dd$n
  
  
  ## constrained parameters
  xfix <- setdiff(colnames(data)[(colnames(data) %in% lava::parlabels(x, exo = TRUE))], lava::latent(x))
  fix <- ifelse(length(xfix) > 0, FALSE, TRUE)
  
  ## reindex the model
  x <- lava::fixsome(x, measurement.fix = fix, S = S, mu = mu, 
                     n = n, debug = FALSE)
  
  if (length(xfix) > 0) {
    index(x) <- lava::reindex(x, sparse = optim$sparse, zeroones = TRUE, 
                              deriv = TRUE)
  }
  else {
    x <- lava::updatelvm(x, sparse = optim$sparse, zeroones = TRUE, 
                         deriv = TRUE, mean = TRUE)
  }
  
  ## number of parameters
  nparall <- lava::index(x)$npar + ifelse(optim$meanstructure, lava::index(x)$npar.mean + lava::index(x)$npar.ex, 0)
  
  myparnames <- coef(x, mean = TRUE)
  
  ## potential input
  paragree <- FALSE
  paragree.2 <- c()
  if (!is.null(optim$start)) {
    paragree <- myparnames %in% names(optim$start)
    paragree.2 <- names(optim$start) %in% myparnames
  }
  
  ## start function
  if (is.null(optim$starterfun) && lava.options()$param != "relative"){optim$starterfun <- startvalues0}
  
  ##
  start <- suppressWarnings(do.call(optim$starterfun, 
                                    list(x = x, S = S, mu = mu, debug = lava.options()$debug, 
                                         silent = TRUE, data = data, ...)))
  
  ## update parameters
  
  if (!is.null(x$expar) && length(start) < nparall) {
    ii <- which(index(x)$e1 == 1)
    start <- c(start, structure(unlist(x$expar[ii]), 
                                names = names(x$expar)[ii]))
  }
  if (length(paragree.2) > 0) {
    start[which(paragree)] <- optim$start[which(paragree.2)]
  }
  
  ## export
  return(start)
  
}





