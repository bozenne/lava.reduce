#' #' @title Extract moments from a reduced lvm
#' #' @description Extract moments from a reduced lvm
#' #' 
#' #' @param x latent variable model
#' #' @param p Parameter vector used to calculate statistics
#' #' @param ... additional arguments to be passed to the low level functions
#' #' 
#' 
#' #' @export
#' coef.lvmfit.reduced <- function (object, ...) {
#'   
#'   object.coef <- callS3methodParent(object, FUN = "coef", class = "lvmfit.reduced")
#'   object$xfull <- reduce2lvm(Model(object))
#'   
#'   
#'   object$xfull
#'   res <- callS3methodParent(object, FUN = "coef", class = "lvmfit.reduced", ...)
#'   
#'   
#'   
#'   stdCoef  <- lava:::stdcoef(xfull, p = object.coef[coef(xfull)])
#'   
#'   
#'   
#'   # res[,2] <- stdCoef$vstar
#'   
#'   return(res)
#' }