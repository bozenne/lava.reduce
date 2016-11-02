
#### timing
test <- FALSE

if(test == TRUE){
  
  method <- lava:::gaussian_method.lvm
  timout <- 12
  
  system.time( 
    LVM <- try(R.utils:::evalWithTimeout(estimate(m, data = d, control = list(start = start, method = method), estimator = "gaussian2"), timeout = timout))
  )
  # gaussian2_method.lvm
  system.time(
    LVM.red <- try(R.utils:::evalWithTimeout(estimate(m.red, data = d, control = list(start = start[coef(m.red)], method = method), estimator = "gaussian2"), timeout = timout))
  )
  
  
  system.time( 
    pLVM <- try(R.utils:::evalWithTimeout(estimate(pm, data = d, lambda1 = lambda1, control = list(start = start, iter.max = 10, method = method), estimator = "gaussian2"), timeout = timout))
  )
  # gaussian2_method.lvm
  system.time( 
    pLVM.red <- try(R.utils:::evalWithTimeout(estimate(pm.red, data = d, lambda1 = lambda1, control = list(start = start[coef(m.red)], iter.max = 10, method = method), estimator = "gaussian2"), timeout = timout))
  )
  
}

