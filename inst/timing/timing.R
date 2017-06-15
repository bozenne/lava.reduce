library(ggplot2)
library(lava.reduce)
library(rbenchmark)
library(data.table)
library(R.utils)

grid <- expand.grid(n = c(50, 500, 1000),
                    p = c(5,10,25,50,100,200))
n.grid <- NROW(grid)
dt <- NULL
n.rep <- 1

# iter.grid <- 1
time.out <- 15

set.seed(10)
pb <- txtProgressBar(max = n.grid)
for(iter.grid in 1:n.grid){
  
  if(grid[iter.grid,"p"]>=grid[iter.grid,"n"]){next}
  
  # model
  m <- lvm()
  m <- regression(m,y='y1',x='x'%++%1:grid[iter.grid,"p"])
  mR <- reduce(m)
  
  # simul
  d <- sim(m,grid[iter.grid,"n"])
  

  tps1 <- system.time(
    e1 <- tryCatch(evalWithTimeout(estimate(m,d), timeout = time.out), silent = TRUE, 
                   error = function(cond){return(list(opt = list(convergence=1, iterations = NA), control = list(method = "nlminb2")))})
  )
  tps2 <- system.time(
    e2 <- tryCatch(evalWithTimeout(estimate(m,d, estimator = "gaussian2"), timeout = time.out), silent = TRUE, 
                   error = function(cond){return(list(opt = list(convergence=1, iterations = NA), control = list(method = "NR")))})
  )
  tps3 <- system.time(
    e3 <- tryCatch(evalWithTimeout(estimate(mR,d), timeout = time.out), silent = TRUE, 
                   error = function(cond){return(list(opt = list(convergence=1, iterations = NA), control = list(method = "nlminb2")))})
  )
  tps4 <- system.time(
    e4 <- tryCatch(evalWithTimeout(estimate(mR,d, estimator = "gaussian2", control = list(method = "NR")), timeout = time.out), silent = TRUE, 
                   error = function(cond){return(list(opt = list(convergence=1, iterations = NA), control = list(method = "NR")))})
  )
  if(is.null(e2$opt$convergence)){
    e2$opt$convergence <- if(e2$opt$iterations<lava.options()$iter.max){0}else{1}
  }
  if(is.null(e4$opt$convergence)){
    e4$opt$convergence <- if(e4$opt$iterations<lava.options()$iter.max){0}else{1}
  }
    
  dt <- rbind(dt,
              data.table(n = grid[iter.grid,"n"], p = grid[iter.grid,"p"],
                         method = c("estimator=gaussian","estimator=gaussian2","estimator=gaussian2LP","estimator=gaussian2LP-NR"), 
                         optim = c(e1$control$method,e2$control$method, e3$control$method, e4$control$method),
                         elapsed = c(tps1[3], tps2[3], tps3[3], tps4[3]), 
                         iterations = c(e1$opt$iterations,e2$opt$iterations, e3$opt$iterations, e4$opt$iterations),
                         cv = c(e1$opt$convergence==0,e2$opt$convergence==0, e3$opt$convergence==0, e4$opt$convergence==0))
  )
  
  setTxtProgressBar(pb, iter.grid)
  
}

#### plot 

ggTime <- ggplot(dt, aes(x = p, y = elapsed, group = method, color = method))
ggTime <- ggTime + geom_point(aes(shape = cv), size = 2) + geom_line()
ggTime <- ggTime + facet_wrap(~n)
ggTime <- ggTime + geom_abline(intercept = time.out, slope = 0, color = "black", linetype = 2)
ggTime
ggsave(ggTime, file = "lava.reduce/inst/timing/plotTime.png")
ggsave(ggTime, file = "lava.reduce/inst/timing/plotTime.svg")

ggIter <- ggplot(dt, aes(x = p, y = iterations, group = method, color = method))
ggIter <- ggIter + geom_jitter(aes(shape = cv), size = 2, width = 0, height = 0) + geom_line()
ggIter <- ggIter + facet_wrap(~n)
ggIter
ggsave(ggIter, file = "lava.reduce/inst/timing/plotIter.png")
ggsave(ggIter, file = "lava.reduce/inst/timing/plotIter.svg")



# {{{ Example

if(FALSE){
    
    library(lavaReduce)

    n <- 100
    p <- 10

    m <- lvm()
    m <- regression(m,y='y1',x='x'%++%1:p)
    mR <- reduce(m, clean = TRUE, rm.exo = TRUE)
    
    mR2 <- lvm.reduced()  # FASTER
    mR2 <- regression(mR2,y='y1',x='x'%++%1:p, reduce = TRUE)
    expect_equal(mR,mR2)

    mR$M
    mR2$M

    # simul
    d <- sim(m,n)
    d[,lp(mR)] <- 0

    #### gradient
    e1 <- estimate(m,d)
    start1 <- coef(e1)[coef(mR)]
    start2 <- coef(e1)

}


# }}}
