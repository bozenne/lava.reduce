context("#### reduce2lvm #### \n")

m <- lvm(c(Y1~X1+X2+X3+eta,Y2~eta,Y3~Z1+X1+eta))
latent(m) <- ~eta
mR <- reduce(m)
mIR <- reduce2lvm(mR)

test_that("reduce2lvm.lvm", {
  expect_equal(mIR$M[rownames(m$M),colnames(m$M)],m$M)
  expect_equal(mIR$par[rownames(m$par),colnames(m$par)],m$par)
  expect_equal(mIR$cov[rownames(m$cov),colnames(m$cov)],m$cov)
  expect_equal(mIR$covpar[rownames(m$covpar),colnames(m$covpar)],m$covpar)
  expect_equal(mIR$fix[rownames(m$fix),colnames(m$fix)],m$fix)
  expect_equal(mIR$covfix[rownames(m$covfix),colnames(m$covfix)],m$covfix)
  expect_equal(mIR$latent,m$latent)
  expect_equal(mIR$mean[names(m$mean)],m$mean)
  # expect_equal(mIR$index,m$index)
  expect_equal(mIR$exogenous,m$exogenous)
  expect_equal(mIR$constrain,m$constrain)
  expect_equal(mIR$constrainY,m$constrainY)
  # expect_equal(mIR$attributes,m$attributes)
})

m <- lvm(c(Y1~X1+X2+X3,Y2~X2+Z3))
mR <- reduce(m)
set.seed(10)
d <- sim(m, 1e3)
e <- estimate(m, d)
eR <- estimate(mR, d, control = list(constrain = TRUE))


test_that("reduce2lvm.lvm", {
  mRfit <- reduce2lvm(eR)
  mfit <- e$model
  expect_equal(mRfit$M[rownames(mfit$M),colnames(mfit$M)],mfit$M)
  expect_equal(mRfit$par[rownames(mfit$par),colnames(mfit$par)],mfit$par)
  expect_equal(mRfit$cov[rownames(mfit$cov),colnames(mfit$cov)],mfit$cov)
  expect_equal(mRfit$covpar[rownames(mfit$covpar),colnames(mfit$covpar)],mfit$covpar)
  expect_equal(mRfit$fix[rownames(mfit$fix),colnames(mfit$fix)],mfit$fix)
  expect_equal(mRfit$covfix[rownames(mfit$covfix),colnames(mfit$covfix)],mfit$covfix)
  expect_equal(mRfit$latent,mfit$latent)
  expect_equal(mRfit$mean[names(mfit$mean)],mfit$mean)
  # expect_equal(mIR$index,mfit$index)
  expect_equal( mfit$index$P[rownames(mfit$index$P), colnames(mfit$index$P)],mfit$index$P)
  expect_equal(mRfit$exogenous,mfit$exogenous)
  expect_equal(mRfit$constrain,mfit$constrain)
  expect_equal(mRfit$constrainY,mfit$constrainY)
})

test_that("moments", {
  momR <- moments(eR)
  mom <- moments(e)
  expect_equal(momR$Cfull[rownames(mom$Cfull),colnames(mom$Cfull)],mom$Cfull, tol = 1e-5)  
})

test_that("iid", {
  iidR <- iid(eR)
  attr(iidR,"bread") <- NULL
  iid0 <- iid(e)
  attr(iid0,"bread") <- NULL
  expect_equal(iidR[,colnames(iid0)], iid0, tol = 1e-6)  
})

