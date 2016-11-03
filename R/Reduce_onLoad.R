'.onLoad' <- function(lib, pkg="lava.reduce") {
  
  lava::addhook("lava.reduce.estimate.hook", hook = "estimate.hooks")
  lava::addhook("lava.reduce.post.hook", hook = "post.hooks")
 
}

'.onAttach' <- function(lib, pkg="lava.reduce") {
  desc <- utils::packageDescription(pkg)
  packageStartupMessage(desc$Package, " version ",desc$Version)
}


deriv.lvm <- get("deriv.lvm", envir = asNamespace("lava"), inherits = FALSE)

estimate.lvm <- get("estimate.lvm", envir = asNamespace("lava"), inherits = FALSE)

`latent<-.lvm` <- get("latent<-.lvm", envir = asNamespace("lava"), inherits = FALSE)
  
regression.lvm <- get("regression.lvm", envir = asNamespace("lava"), inherits = FALSE)
  


