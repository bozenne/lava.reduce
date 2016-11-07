# library(testthat)
# library(lava.reduce)
# library(lava)

context("#### destructors #### \n")

#### regresssion ####
m <- lvm.reduced()
m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
m <- regression(m, x="z",y="y")

## cancel
test_that("cancel", {
  mc <- cancel(m, coef(m)[2])
  expect_false(coef(m)[2] %in% coef(mc))
  expect_true(all(coef(m)[-2] %in% coef(mc)))
  
  mc <- cancel(m, coef(m)[4])
  expect_false(coef(m)[4] %in% coef(mc))
  
  mm <- cancel(m, coef(m)[4], restaure = TRUE)
  expect_true(coef(m)[4] %in% coef(mm))
  expect_false(coef(m)[4] %in% lp(mm, type = "link"))
  
  mm <- cancel(mm, coef(m)[4:13])
  expect_false(any(coef(m)[4:13] %in% coef(mm)))
  
  mm <- cancel(m, coef(m)[4:13], restaure = TRUE)
  expect_true(all(coef(m)[4:13] %in% coef(mm)))
  expect_false("lvm.reduced" %in% class(mm))
})

## kill
m <- lvm.reduced()
m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)

test_that("kill - regression",{
  m1 <- m
  kill(m1) <- ~LPy
  expect_false("lvm.reduced" %in% class(m1))
  
  m2 <- m
  kill(m2) <- ~x1
  expect_true("lvm.reduced" %in% class(m2))
  expect_false(coef(m)[3] %in% coef(m2))
  
})

test_that("kill - LVM",{
  m <- lvm.reduced(c(Y1,Y2,Y3)~eta+Z)
  m <- regression(m, x="x1",y="Y3")
  m <- regression(m, x=paste0("x",1:10),y="Y1", reduce = TRUE)
  m <- regression(m, x=paste0("x",1:10),y="Y2", reduce = TRUE)
  
  m1 <- m
  kill(m1) <- ~LPY1+LPY2
  expect_false("lvm.reduced" %in% class(m1))
  
  m2 <- m
  kill(m2) <- ~LPY1
  expect_equal(lp(m2), lp(m, lp = 2))
  
  m3 <- m
  kill(m3) <- ~Z
  expect_false("Z" %in% vars(m3))
  expect_true("lvm.reduced" %in% class(m3))
  
  m4 <- m
  kill(m4) <- ~x1
  expect_false("x1" %in% vars(m4))
  expect_false("x1" %in% lp(m4, type = "x"))
  expect_true("lvm.reduced" %in% class(m4))
  
  m5 <- m
  kill(m5) <- ~x5
  expect_false("x5" %in% vars(m5))
  expect_false("x5" %in% lp(m5, type = "x"))
  expect_true("lvm.reduced" %in% class(m5))
  
})