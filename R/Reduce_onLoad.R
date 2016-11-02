'.onLoad' <- function(lib, pkg="lava.reduce") {
  
  lava::addhook("lava.reduce.estimate.hook", hook = "estimate.hooks")
  lava::addhook("lava.reduce.post.hook", hook = "post.hooks")
  
}

'.onAttach' <- function(lib, pkg="lava.reduce") {
  desc <- utils::packageDescription(pkg)
  packageStartupMessage(desc$Package, " version ",desc$Version)
}