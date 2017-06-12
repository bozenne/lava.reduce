lavaReduce.env <- new.env()
assign("clean.hooks",c(),envir=lavaReduce.env)
assign("reduce.hooks",c(),envir=lavaReduce.env)

'.onLoad' <- function(lib, pkg="lavaReduce") {
 
    lava::addhook("lavaReduce.estimate.hook", hook = "estimate.hooks")
    lava::addhook("lavaReduce.post.hook", hook = "post.hooks")
    lava::addhook("lavaReduce.remove.hook", hook = "remove.hooks")
    lava::addhook("lavaReduce.cancel.hook", hook = "cancel.hooks")
     
    lava::lava.options(estimator.default.reduce = "2",  # switch from gaussian moment to gaussian2 moment
                       init.restaure.reduce = TRUE) # initialise using all variables (else set the variables of the linear predictors to 0)

    addhook_lavaReduce("lavaReduce.clean.hook", hook = "clean.hooks")
}


'.onAttach' <- function(lib, pkg="lavaReduce") {
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
#' @param ... for compatibility with lava

#' @rdname hook.reduce
#' @export
gethook_lavaReduce<- function (hook, ...){
    get(hook, envir = lavaReduce.env)
}

#' @rdname hook.reduce
#' @export
addhook_lavaReduce <- function (x, hook, ...){
    newhooks <- unique(c(gethook_lavaReduce(hook), x))
    assign(hook, newhooks, envir = lavaReduce.env)
    invisible(newhooks)
}
