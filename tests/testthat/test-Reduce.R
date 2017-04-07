# library(testthat)
# library(lava.reduce)
# butils.base::package.source("lava.reduce", Rcode = TRUE, RorderDescription = FALSE)

gaussian1LP_gradient.lvm <- lava.reduce:::gaussian1LP_gradient.lvm
  
context("#### Reduce #### \n")
lava.options(symbols = c(";","<->"))

m <- lvm()
regression(m) <- y1 ~ x1 + x2
m.red <- reduce(m)

test_that("reduce twice",{
  expect_equal( reduce(m.red) , m.red )
})

iter.max <- 1

context("#### Regression #### \n")

m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:2)
m <- regression(m,y='y1',x='z'%++%1:5)
m.red <- reduce(m)

## simul
set.seed(10)
d <- sim(m,5e2)

suppressWarnings(
  start <- coef(estimate(m, data = d, control = list(iter.max = 0), quick = FALSE))
)

## models

## tests moment
test_that("Regression: moment reduce", {
  e <- estimate(m, d, estimator = "gaussian1")
  index <- match(coef(m.red),coef(m))
  
  dLP <- d
  dLP[,lp(m.red)] <- 0
  
  g1 <- gaussian1LP_gradient.lvm(m.red, p = start[index], data = dLP, indiv = FALSE)
  g2 <- lava:::gaussian1_gradient.lvm(x = m, data=d, p=start, S = e$S, n = e$data$n, mu = e$mu)
  expect_equal(unname(g1), g2[index])

  # gaussianLP_logLik.lvm <- lava.reduce:::gaussianLP_logLik.lvm
  # gaussianLP_gradient.lvm <- lava.reduce:::gaussianLP_gradient.lvm
  # gaussianLP_score.lvm <- lava.reduce:::gaussianLP_score.lvm
  # gaussianLP_hessian.lvm <- lava.reduce:::gaussianLP_hessian.lvm
  # gaussian1LP_hessian.lvm <- lava.reduce:::gaussian1LP_hessian.lvm
  # gaussian2LP_hessian.lvm <- lava.reduce:::gaussian2LP_hessian.lvm
  # 
  # checkMoment(m.red, e) # the data must be converted to matrix + add LP
    
})

## tests estimation
test_that("Regression: lvm vs lvm.reduce", {
  
  method <- lava:::gaussian_method.lvm
  suppressWarnings(
    LVM1 <- estimate(m, data = d,  control = list(iter.max = iter.max, start = start, method = method), estimator = "gaussian1")
  )
  suppressWarnings(
    LVM1.red <- estimate(m.red, data = d, control = list(iter.max = iter.max, start = start[coef(m.red)], method = method), estimator = "gaussian1")
  )
  
  expect_equal(coef(LVM1),coef(LVM1.red)[names(coef(LVM1))])
  
  method <- lava:::gaussian2_method.lvm
  suppressWarnings( 
    LVM2 <- estimate(m, data = d, control = list(iter.max = iter.max, start = start, method = method), estimator = "gaussian2")
  )
  suppressWarnings(
    LVM2.red <- estimate(m.red, data = d, control = list(iter.max = iter.max, start = start[coef(m.red)], method = method), estimator = "gaussian2")
  )
  expect_equal(coef(LVM2),coef(LVM2.red)[names(coef(LVM2))])
})

context("#### Latent variable model #### \n")
# m <- lvm()
# m <- regression(m,y=c('y1','y2','y3','y4'),x='eta')
# m <- regression(m,y=c('y2','y3'),x='x'%++%1:5)
# latent(m) <- ~eta
# m <- regression(m,y=c('y4','y2'),x='z'%++%1:2)
# covariance(m) <- y2~y1
m <- lvm()
m <- regression(m,y=c('y1','y2','y3'),x='eta')
m <- regression(m,y=c('y2'),x='x'%++%1)
latent(m) <- ~eta
m <- regression(m,y=c('y3'),x='z'%++%1)


m.red1 <- reduce(m)
m.red2 <- reduce(m, endo = "y3")

## simul
set.seed(10)
d <- sim(m,5e2, latent = FALSE)
start <- setNames(rep(0, length(coef(m))), coef(m))
start[grep(lava.options()$symbol[1], names(start))] <- 1
suppressWarnings(
  startLVM <- coef(estimate(m, data = d, control = list(iter.max = 0)))
)
start[names(startLVM)] <- startLVM
# names(start) <- gsub("~~",",",names(start))
# names(startLVM) <- gsub("~~",",",names(startLVM))

e <- estimate(m, d, estimator = "gaussian1")
index1 <- match(coef(m.red1),coef(m))
index2 <- match(coef(m.red2),coef(m))

## tests moment
test_that("LVM: moment reduce", {
 
  g1 <- gaussian1LP_gradient.lvm(m.red1, p = start[index1], data = d, indiv = FALSE)
  g2 <- gaussian1LP_gradient.lvm(m.red2, p = start[index2], data = d, indiv = FALSE)
  g3 <- lava:::gaussian1_gradient.lvm(x = e$model, data=d, p=startLVM, S = e$S, n = e$data$n, mu = e$mu)
  g3 <- setNames(g3, names(startLVM))
  
  expect_equal(g1[intersect(names(g1),names(g3))], 
               g3[intersect(names(g1),names(g3))])
  expect_equal(g2[intersect(names(g2),names(g3))], 
               g3[intersect(names(g2),names(g3))])
})

## tests estimation
method <- lava:::gaussian_method.lvm
suppressWarnings(
  start1 <- estimate(m.red1, data = d, control = list(iter.max = 0), quick = TRUE)
)
start1 <- start[names(start1)]

suppressWarnings(
  start2 <- estimate(m.red2, data = d, control = list(iter.max = 0), quick = TRUE)
)
start2 <- start[names(start2)]

test_that("LVM: lvm vs lvm.reduce (gaussian 1)", {
  suppressWarnings(
    LVM1.red <- estimate(m.red1, data = d, 
                          control = list(iter.max = iter.max, start = start1, method = method), quick = TRUE, 
                          estimator = "gaussian1")
  )
  
 
  suppressWarnings(
    LVM2.red <- estimate(m.red2, data = d, 
                         control = list(iter.max = iter.max, start = start2, method = method), quick = TRUE, 
                         estimator = "gaussian1")
  )
  
  suppressWarnings(
    LVM3 <- estimate(m, data = d,  
                     control = list(iter.max = iter.max, start = start, method = method),
                     estimator = "gaussian1")
  )
  
  expect_equal(LVM1.red[intersect(names(LVM1.red),names(coef(LVM3)))],
               coef(LVM3)[intersect(names(LVM1.red),names(coef(LVM3)))]
               )
  
  expect_equal(LVM2.red[intersect(names(LVM1.red),names(coef(LVM3)))],
               coef(LVM3)[intersect(names(LVM1.red),names(coef(LVM3)))]
  )
})

test_that("LVM: lvm vs lvm.reduce (gaussian 2)", {
  
  suppressWarnings(
    LVM1.red <- estimate(m.red1, data = d, 
                         control = list(iter.max = iter.max, start = start1, method = method), quick = TRUE, 
                         estimator = "gaussian2")
  )
  
  
  suppressWarnings(
    LVM2.red <- estimate(m.red2, data = d, 
                         control = list(iter.max = iter.max, start = start2, method = method), quick = TRUE, 
                         estimator = "gaussian2")
  )
  
  suppressWarnings(
    LVM3 <- estimate(m, data = d,  
                     control = list(iter.max = iter.max, start = start, method = method),
                     estimator = "gaussian2")
  )
  
  expect_equal(LVM1.red[intersect(names(LVM1.red),names(coef(LVM3)))],
               coef(LVM3)[intersect(names(LVM1.red),names(coef(LVM3)))]
  )
  
  expect_equal(LVM2.red[intersect(names(LVM1.red),names(coef(LVM3)))],
               coef(LVM3)[intersect(names(LVM1.red),names(coef(LVM3)))]
  )
})


#### case where the linear predictor is first 
m <- lvm()
m <- regression(m,y=c('y1','y2','y3'),x='eta')
m <- regression(m,y=c('y1'),x='x'%++%1)
latent(m) <- ~eta
m <- regression(m,y=c('y2'),x='z'%++%1)

## simul
set.seed(10)
d <- sim(m,5e2, latent = FALSE)

m <- reduce(m)

suppressWarnings(
LVM1.red <- estimate(m, data = d, 
                     control = list(iter.max = 1), quick = TRUE, 
                     estimator = "gaussian1")
)
# automatic switch to y3 as a reference


#### can reduce with only one variable
m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:1)
test_that("reduce only one variable", {
  mR <- reduce(m)
})


####  to fix ####

if(FALSE){
set.seed(10)

# model
m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:5)
mR <- reduce(m)

# simul
d <- sim(m,50)

e <- estimate(mR,d, estimator = "gaussian2", control = list(trace = 2))

estimate(mR,d, estimator = "gaussian2", control = list(trace = 2, method = "nlminb2"))
estimate(m,d, estimator = "gaussian2", control = list(trace = 2, method = "nlminb2"))

estimate(mR,d, estimator = "gaussian2", control = list(trace = 2, method = "NR"))
estimate(m,d, estimator = "gaussian2", control = list(trace = 2, method = "NR"))

estimate(mR,d, estimator = "gaussian2", control = list(trace = 2, method = "nlminb1"))
estimate(m,d, estimator = "gaussian2", control = list(trace = 2, method = "nlminb1"))
}
