library(ggplot2)
library(lavaReduce)
library(rbenchmark)
library(data.table)
library(R.utils)

n <- 1000
seq_p <- c(5,10,15,20,25,40,50,75,100)
n.p <- length(seq_p)

dt <- NULL
pb <- txtProgressBar(max = n.p)

for(iter_p in 1:n.p){
  # model
  m <- lvm()
  m <- regression(m,y='y1',x='x'%++%1:seq_p[iter_p])
  mR <- reduce(m, clean = TRUE)
  
  # simul
  d <- sim(m,n)
  d[,lp(mR)] <- 0
  #### gradient
  e1 <- estimate(m,d)
  start1 <- coef(e1)[coef(mR)]
  start2 <- coef(e1)
  
  gaussianLP_score.lvm(x = mR, data = as.matrix(d), p = start1, method = "R")
  gaussianLP_score.lvm(x = mR, data = as.matrix(d), p = start1, implementation = "cpp")
  ben <- rbenchmark::benchmark(reduceCpp = ,
                               reduceR = ,
                               lava = lava:::gaussian_gradient.lvm(x = e1$model, data = d, p = start2, n = e1$data$n, mu = e1$mu, S = e1$S),  # default
                               replications = 10
  )
  gaussianLP_score.lvm
  
  
  
  # FUN <- function(){gaussianLP_gradient.lvm(x = mR, data = as.matrix(d), p = start1, method = "cpp")}
  # butils::profileCode(FUN)
  ben <- rbenchmark::benchmark(reduceCpp = gaussianLP_gradient.lvm(x = mR, data = as.matrix(d), p = start1, method = "cpp"),
                               reduceR = gaussianLP_gradient.lvm(x = mR, data = as.matrix(d), p = start1, method = "R"),
                               lava = lava:::gaussian_gradient.lvm(x = e1$model, data = d, p = start2, n = e1$data$n, mu = e1$mu, S = e1$S),  # default
                               replications = 10
  )
  
  setTxtProgressBar(pb, value = iter_p)  
  
  dt <- rbind(dt,
              cbind(as.data.table(ben[c("test","elapsed")]), p = seq_p[iter_p])
  )
}


dt[, test := factor(test, labels = c("lava.reduceCpp","lava.reduceR","lava"))]


gg <- ggplot(dt, aes(x = p, y = elapsed, color = test, group = test))
gg <- gg + geom_line() + geom_point() + scale_color_hue("gradient")
gg


ggsave("inst/timing/plotTimeGrad.svg")
ggsave("inst/timing/plotTimeGrad.png")

system.time(
  G1 <- gaussian1LP_gradient.lvm(x = mR, data = d, p = coef(e1))
)
system.time(
  G2 <- lava:::gaussian1_gradient.lvm(x = e1$model, data = d, p = coef(e1), n = e1$data$n, mu = e1$mu, S = e1$S)
)

system.time(
  H1 <- gaussian1LP_hessian.lvm(x = mR, data = d, p = coef(e1))
)
system.time(
  H2 <- lava:::gaussian1_hessian.lvm(x = e1$model, data = d, p = coef(e1), n = e1$data$n, mu = e1$mu, S = e1$S)
)


