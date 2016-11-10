'.onLoad' <- function(lib, pkg="lava.reduce") {
  
  lava::addhook("lava.reduce.estimate.hook", hook = "estimate.hooks")
  lava::addhook("lava.reduce.post.hook", hook = "post.hooks")
 
  lava::lava.options(estimator.default.reduce = "2",  # switch from gaussian moment to gaussian2 moment
                     init.restaure.reduce = TRUE) # initialise using all variables (else set the variables of the linear predictors to 0)
}

'.onAttach' <- function(lib, pkg="lava.reduce") {
  desc <- utils::packageDescription(pkg)
  packageStartupMessage(desc$Package, " version ",desc$Version)
}

procdata.lvm <- get("procdata.lvm", envir = asNamespace("lava"), inherits = FALSE)
