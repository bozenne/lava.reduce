### Lava_cleanExo.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar 10 2017 (10:22) 
## Version: 
## last-updated: apr  4 2017 (09:28) 
##           By: Brice Ozenne
##     Update #: 87
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ doc
#' @title Simplify a lvm object
#' @name clean
#' @description Remove variables with no link and simplify the class of the lvm object
#' 
#' @param x \code{lvm}-object
#' @param rm.exo  should the exogenous variables with no links be removed from the object ? 
#' @param rm.endo  should the endogenous variables with no links be removed from the object ?
#' @param rm.latent  should the latent variables with no links be removed from the object ?
#' @param rm.lp  should the linear predictors with no links be removed from the object ?
#' @param simplify.reduce  should the class \code{lvm.reduced} be droped if there is no linear predictor in the object?
#' @param simplify  should the class of the object be simplified ? Overwrite the simplify.x arguments.
#' @param ... additional arguments to lower level functions
#'
#' @details
#' simplify means remove the class \code{"lavaReduce"} if no linear predictor is in the object.
#' 
#' @examples 
#' 
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:5),y="y1")
#' m <- regression(m, x=paste0("x",1:5),y="y2")
#' covariance(m) <- y1~y2
#'
#' cancel(m) <- y1 ~ x1
#' cancel(m) <- y2 ~ x1
#' clean(m)
#'
#' m <- lvm(y1 ~ eta + x1, y2 ~ eta, y3 ~ eta + x2)
#' latent(m) <- ~eta
#' clean(m)
#' m
#' cancel(m) <- y1 ~ eta
#' cancel(m) <- y2 ~ eta
#' cancel(m) <- y3 ~ eta
#' clean(m)
#' 
#' @export
`clean` <-
  function(x,...) UseMethod("clean")
# }}}

# {{{ clean.lvm
#' @rdname clean
#' @export
clean.lvm <- function(x, rm.exo = TRUE, rm.endo = TRUE, rm.latent = TRUE, ...){

    myhooks <- lavaReduce::gethook_lavaReduce("clean.hooks")
    for (f in myhooks) {
        res <- do.call(f, list(x=x, rm.exo=rm.exo,...))
        if("x" %in% names(res)){ x <- res$x }
        if("rm.exo" %in% names(res)){ rm.exo <- res$rm.exo }
        if("rm.endo" %in% names(res)){ rm.endo <- res$rm.endo }
    }  
    
    var.latent <- latent(x)
    var.exogenous <- exogenous(x)
    var.endogenous <- endogenous(x)

    M.reg <- x$index$A
    M.cov <- x$index$P - diag(diag(x$index$P))
    
    varClean <- NULL
    if(rm.exo && length(var.exogenous) > 0){
        indexClean.reg <- which(rowSums(M.reg[var.exogenous,,drop = FALSE]!=0)==0)
        indexClean.cov <- which(rowSums(M.cov[var.exogenous,,drop = FALSE]!=0)==0)
        indexClean <- intersect(indexClean.reg, indexClean.cov)
        varClean <- c(varClean, var.exogenous[indexClean])
    }
    if(rm.endo && length(var.endogenous)>0){
        indexClean.reg <- which(colSums(M.reg[,var.endogenous,drop = FALSE]!=0)==0)
        indexClean.cov <- which(colSums(M.cov[,var.endogenous,drop = FALSE]!=0)==0)
        indexClean <- intersect(indexClean.reg, indexClean.cov)
        varClean <- c(varClean, var.endogenous[indexClean])
    }
    if(rm.latent && length(var.latent)>0){
        indexClean.Rreg <- which(rowSums(M.reg[var.latent,,drop = FALSE]!=0)==0)
        indexClean.Rcov <- which(rowSums(M.cov[var.latent,,drop = FALSE]!=0)==0)
        indexClean.Creg <- which(colSums(M.reg[,var.latent,drop = FALSE]!=0)==0)
        indexClean.Ccov <- which(colSums(M.cov[,var.latent,drop = FALSE]!=0)==0)
        indexClean <- intersect(intersect(indexClean.Rreg, indexClean.Rcov),
                                intersect(indexClean.Creg, indexClean.Ccov))
        varClean <- c(varClean, var.latent[indexClean])
    }
    
    if(length(varClean)>0){
        x <- kill(x, value =  varClean, ...)
    }
    return(x)
}
# }}}

# {{{ lavaReduce.clean.hook
#' @rdname clean
#' @export
lavaReduce.clean.hook <- function(x, rm.exo,
                                   rm.lp = TRUE, simplify.reduce = TRUE, simplify, ...){

    if("lvm.reduced" %in% class(x)){
        if(!missing(simplify)){
            simplify.reduce <- simplify
        }

        ## remove empty lp
        n.link <- lp(x, type = "n.link")
        if(rm.lp && any(n.link==0)){ 
            name.lp <- names(n.link)[which(n.link==0)]
            x <- kill(x, value = name.lp)
            x$lp[lp(x, lp = name.lp, type = "endo")] <- NULL
        }

        ## remove useless exogeneous variables
        if(length(lp(x, type = "link", format = "vector")) == 0 ){
            ## lvm.reduced object with no linear predictor
            if(simplify.reduce){
                class(x) <- setdiff(class(x), "lvm.reduced")
            }
        }else{
            ## lvm.reduced object with linear predictors
            var.exogenous <- exogenous(x)
            M.reg <- x$index$A
            M.cov <- x$index$P - diag(diag(x$index$P))

            if(rm.exo && length(var.exogenous) > 0){
                indexClean.reg <- which(rowSums(M.reg[var.exogenous,,drop = FALSE]!=0)==0)
                indexClean.cov <- which(rowSums(M.cov[var.exogenous,,drop = FALSE]!=0)==0)
                # indexClean.lp <- which(var.exogenous %in% lp(x, type = "x") == FALSE)
                indexClean <- intersect(indexClean.reg, intersect(indexClean.cov, indexClean.lp))
                if(length(indexClean)>0){
                    x <- kill(x, value =  var.exogenous[indexClean])
                }
            }

            rm.exo <- FALSE
         
        }
    }

    return(list(x=x,
                rm.exo=rm.exo))

}
# }}}

#----------------------------------------------------------------------
### Lava_cleanExo.R ends here
