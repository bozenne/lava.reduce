lava.reduce.env <- new.env()
assign("clean.hooks",c(),envir=lava.reduce.env)
assign("reduce.hooks",c(),envir=lava.reduce.env)

'.onLoad' <- function(lib, pkg="lava.reduce") {
 
    lava::addhook("lava.reduce.estimate.hook", hook = "estimate.hooks")
    lava::addhook("lava.reduce.post.hook", hook = "post.hooks")
    lava::addhook("lava.reduce.remove.hook", hook = "remove.hooks")
    lava::addhook("lava.reduce.cancel.hook", hook = "cancel.hooks")
     
    lava::lava.options(estimator.default.reduce = "2",  # switch from gaussian moment to gaussian2 moment
                       init.restaure.reduce = TRUE) # initialise using all variables (else set the variables of the linear predictors to 0)

    addhook_lava.reduce("lava.reduce.clean.hook", hook = "clean.hooks")
}


'.onAttach' <- function(lib, pkg="lava.reduce") {
  desc <- utils::packageDescription(pkg)
  packageStartupMessage(desc$Package, " version ",desc$Version)
}

procdata.lvm <- get("procdata.lvm", envir = asNamespace("lava"), inherits = FALSE)


#' @title Hooks for lava reduce
#' @description Get and add hook for lava reduce
#' @name hook.reduce
#'
#' @param x the function to add to the hook
#' @param hook the name of the hook
#' 

#' @rdname hook.reduce
#' @export
gethook_lava.reduce<- function (hook, ...){
    get(hook, envir = lava.reduce.env)
}

#' @rdname hook.reduce
#' @export
addhook_lava.reduce <- function (x, hook, ...){
    newhooks <- unique(c(gethook_lava.reduce(hook), x))
    assign(hook, newhooks, envir = lava.reduce.env)
    invisible(newhooks)
}
