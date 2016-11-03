# library(testthat)
# library(lava.reduce)
# library(lava)

context("#### initCoef - regression #### \n")

## sim
endo.name <- "y"
exo.name <- paste0("x",1:10)

m <- lvm()
m <- regression(m, x=exo.name,y=endo.name)
mR <- reduce(m)

set.seed(10)
d <- sim(m,1e2)

## test
test_that("initCoef - regression", {
  
  suppressWarnings(
    start1 <- estimate(m,d, control = list(iter.max = 0), quick = TRUE)
  )
  suppressWarnings(
    start2 <- estimate(mR,d, control = list(iter.max = 0), quick = TRUE)
  )
  
  expect_equal(start1, start2[names(start1)])
  
})


context("#### initCoef - LVM #### \n")

endo.name <- list("y1","y2","y3")
exo.name <- list(paste0("x",1:10),c("z",paste0("x",1:3)),paste0("x",1))

m <- lvm()
m <- regression(m, x=c("eta",exo.name[[1]]), y=endo.name[[1]])
m <- regression(m, x=c("eta",exo.name[[2]]), y=endo.name[[2]])
m <- regression(m, x=c("eta",exo.name[[3]]), y=endo.name[[3]])
latent(m) <- "eta"
mR <- reduce(m, endo = c("y2","y3"))

set.seed(10)
d <- sim(m,1e2, latent = FALSE)

## test
test_that("initCoef - regression", {
  
  suppressWarnings(
    start1 <- coef(estimate(m,d, control = list(iter.max = 0), quick = FALSE))
  )
  suppressWarnings(
    start2 <- estimate(mR,d, control = list(iter.max = 0), quick = TRUE)
  )
  
  expect_equal(start1, start2[names(start1)])
  
})