# library(butils.base)
# package.source("lavaReduce", Ccode = TRUE)

context("#### gradient #### \n")

n <- 100
p <- 10

set.seed(10)
m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:p)
mR <- reduce(m)
mR2 <- lvm.reduced()  # FASTER
mR2 <- regression(mR2,y='y1',x='x'%++%1:p, reduce = TRUE)

test_that("2 ways to reduce lvm",{
  expect_equal(coef(mR),coef(mR2))
})
# simul
d <- sim(m,n)
d[,lp(mR)] <- 0

#### gradient
e1 <- estimate(m,d)
start1 <- coef(e1)[coef(mR)]
start2 <- coef(e1)

test_that("Gradient from lava match those of lavaReduce", {
  res1R <- gaussianLP_gradient.lvm(x = mR, data = as.matrix(d), p = start1, implementation = "R")
  res1Cpp <- gaussianLP_gradient.lvm(x = mR, data = as.matrix(d), p = start1, implementation = "Cpp")
  expect_true(all(abs(res1R)<1e-6))
  expect_true(all(abs(res1Cpp)<1e-6))

  res2R <- gaussianLP_gradient.lvm(x = mR, data = as.matrix(d), p = start1+1, implementation = "R")
  res2Cpp <- gaussianLP_gradient.lvm(x = mR, data = as.matrix(d), p = start1+1, implementation = "Cpp")
  res2Lava <- lava:::gaussian_gradient.lvm(x = e1$model, data = d, p = start2+1, n = e1$data$n, mu = e1$mu, S = e1$S)
  
  expect_equal(res2Lava,as.double(res2Cpp[names(start2)]))
})

test_that("Compare results of the estimation",{
  res1 <- estimate(m,d)
  
  res2 <- estimate(mR,d, control = list(constrain = TRUE))
  
  res3 <- estimate(mR2,d, control = list(constrain = TRUE))
  expect_equal(coef(res1),coef(res2)[names(coef(res1))], tol = 1e-5)
  expect_equal(coef(res1),coef(res3)[names(coef(res1))], tol = 1e-5)
  
  xfull <- reduce2lvm(res2)
  sd1 <- lava:::stdcoef(res1)
  sd2 <- lava:::stdcoef(xfull, p = coef(res2)[coef(xfull)])
  expect_equal(sd1,sd2,tol = 1e-5)
})

## problem
# res2$coef != res1$coef