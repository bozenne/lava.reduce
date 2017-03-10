# library(testthat)
# library(lava.reduce)
# library(lava)

context("#### get Vars #### \n")

endo.name <- "y"
exo.name <- paste0("x",1:10)

m <- lvm.reduced()
m <- regression(m, x=exo.name,y=endo.name, reduce = TRUE)

lp.name <- unname(lp(m))
  
  
test_that("vars", {
  
  expect_equal(vars(m), c(endo.name, lp.name))
  expect_equal(vars(m, lp = FALSE), endo.name)
  expect_equal(vars(m, lp = FALSE, xlp = TRUE), c(endo.name, exo.name))
  
})

test_that("exogenous", {
  
  expect_equal(exogenous(m), lp.name)
  expect_equal(exogenous(m, lp = FALSE), character(0))
  expect_equal(exogenous(m, lp = FALSE, xlp = TRUE), exo.name)
  
})

test_that("endogenous", {
  
  expect_equal(endogenous(m), endo.name)
  expect_equal(endogenous(m, lp = FALSE), endo.name)
  expect_equal(endogenous(m, lp = FALSE, xlp = TRUE), endo.name)
  
})

test_that("endogenous", {
  mGS <- lvm()
  mGS <- regression(mGS, x=exo.name,y=endo.name)
  
  expect_equal(unname(coef(m)[c("m1",paste0("e",1:10),"p1")]), unname(coef(mGS)))
  
})


context("#### lp #### \n")

#### regression
endo.name <- "y"
exo.name <- paste0("x",1:10)
m <- lvm.reduced()
m <- regression(m, x=exo.name,y=endo.name, reduce = TRUE)

test_that("lp", {
  
  expect_equal(unname(lp(m)), "LPy")
  
})

#### LVM
endo.name <- list("y1","y2","y3")
exo.name <- list(paste0("x",1:10),c("z",paste0("x",1:3)),paste0("x",1))

m <- lvm.reduced()
latent(m) <- "eta"
m <- regression(m, y=endo.name[[1]], x="eta", reduce = FALSE)
m <- regression(m, y=endo.name[[2]], x="eta", reduce = FALSE)
m <- regression(m, y=endo.name[[3]], x="eta", reduce = FALSE)
m <- regression(m, x=exo.name[[1]], y=endo.name[[1]], reduce = TRUE)
m <- regression(m, x=exo.name[[2]], y=endo.name[[2]], reduce = TRUE)
m <- regression(m, x=exo.name[[3]], y=endo.name[[3]], reduce = FALSE)

lp.name <- unname(lp(m))

test_that("lp - regression = reduce ", {
  
  expect_equal(unname(lp(m, type = "x", lp = lp.name[1])), exo.name[[1]])
  expect_equal(unname(lp(m, type = "x", lp = lp.name[2])), exo.name[[2]])
  
})

#### LVM2
m <- lvm()
m <- regression(m, x=c("eta",exo.name[[1]]), y=endo.name[[1]])
m <- regression(m, x=c("eta",exo.name[[2]]), y=endo.name[[2]])
m <- regression(m, x=c("eta",exo.name[[3]]), y=endo.name[[3]])
latent(m) <- "eta"
mr <- reduce(m)

lp.name <- unname(lp(mr))

test_that("lp - reduce all model", {
  
  expect_equal(sort(unname(lp(mr, type = "x", lp = lp.name[1]))), sort(exo.name[[1]]))
  expect_equal(sort(unname(lp(mr, type = "x", lp = lp.name[2]))), sort(exo.name[[2]]))
  expect_equal(sort(unname(lp(mr, type = "x", lp = lp.name[3]))), sort(exo.name[[3]]))
  
})
